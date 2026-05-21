import streamlit as st
import pandas as pd
import pyteomics
from pyteomics import mzml
import numpy as np
import os
import gc
from io import BytesIO

# ==========================================
# 1. Interface de Utilizador (UI) - Streamlit
# ==========================================
st.set_page_config(page_title="Análise Metabolómica", layout="wide")

st.title("🔬 Processamento de Dados de Metabolómica")
st.markdown("Carregue a tabela de *targets* e defina os caminhos para as bases de dados e ficheiros mzML.")

# Inputs na interface
uploaded_target = st.file_uploader(
    label="1. Escolha a Tabela de Search (Excel)",
    type=["xlsx", "csv"]
)

folder_mzml = st.text_input(
    "2. Caminho para a pasta dos ficheiros MzML", 
    value="/Volumes/HDD_andre/CIIMAR/Metabolomica/metabolomic_data/MzML_files/pos/"
)

path_adducts = st.text_input(
    "3. Caminho para a Tabela de Adutos", 
    #value="/Volumes/HDD_andre/CIIMAR/Metabolomica/Tabela_Adutos_Limpa.xlsx"
    value="Tabela_Adutos_Limpa.xlsx"
)

path_isotopes = st.text_input(
    "4. Caminho para a Tabela de Isótopos", 
    #value="/Volumes/HDD_andre/CIIMAR/Metabolomica/Tabela_Isotopos_Limpa.xlsx"
    value="Tabela_Isotopos_Limpa.xlsx"
)

# ==========================================
# 2. Lógica de Execução
# ==========================================
if st.button("🚀 Iniciar Análise", type="primary"):
    
    if not uploaded_target:
        st.error("Por favor, carregue a Tabela de Search antes de iniciar.")
    elif not os.path.exists(folder_mzml):
        st.error(f"A pasta MzML não foi encontrada: {folder_mzml}")
    else:
        with st.status("A iniciar análise...", expanded=True) as status:
            try:
                # Carregar bases de dados
                st.write("A carregar bases de dados...")
                df_adducts = pd.read_excel(path_adducts)
                df_adducts = df_adducts[df_adducts['Ion_Mode'] == 'Positive']
                df_isotopos = pd.read_excel(path_isotopes)

                # Extrair valores exatos
                st.write("A extrair diferenças isotópicas...")
                diff_13C = df_isotopos.loc[df_isotopos['Isotopo'] == '13C', 'Dif_Massa_Monoisotopico'].values[0]
                diff_15N = df_isotopos.loc[df_isotopos['Isotopo'] == '15N', 'Dif_Massa_Monoisotopico'].values[0]
                diff_18O = df_isotopos.loc[df_isotopos['Isotopo'] == '18O', 'Dif_Massa_Monoisotopico'].values[0]
                diff_34S = df_isotopos.loc[df_isotopos['Isotopo'] == '34S', 'Dif_Massa_Monoisotopico'].values[0]
                diff_37Cl = df_isotopos.loc[df_isotopos['Isotopo'] == '37Cl', 'Dif_Massa_Monoisotopico'].values[0]
                diff_81Br = df_isotopos.loc[df_isotopos['Isotopo'] == '81Br', 'Dif_Massa_Monoisotopico'].values[0]

                diff_Cl2 = diff_37Cl * 2
                diff_Cl3 = diff_37Cl * 3
                diff_Br2 = diff_81Br * 2
                diff_ClBr = diff_37Cl + diff_81Br

                tolerance_rt_minutes = 0.2
                tolerance_ppm = 5.0

                col_13C = f"Isotope_13C (ref:+{diff_13C:.5f} | ab:1.1%)"
                col_15N = f"Isotope_15N (ref:+{diff_15N:.5f} | ab:0.37%)"
                col_18O = f"Isotope_18O (ref:+{diff_18O:.5f} | ab:0.20%)"
                col_34S = f"Isotope_34S (ref:+{diff_34S:.5f} | ab:4.21%)"

                # Ler tabela de targets carregada
                if uploaded_target.name.endswith('.csv'):
                    df_targets = pd.read_csv(uploaded_target)
                else:
                    df_targets = pd.read_excel(uploaded_target)

                targets_by_sample = df_targets.groupby('Sample')
                final_results = []
                
                # Barra de progresso baseada no número de amostras
                progress_bar = st.progress(0)
                total_samples = len(targets_by_sample)

                # Processamento das amostras
                for idx_sample, (sample_name, sample_info) in enumerate(targets_by_sample):
                    file_name = f"{sample_name}_03.mzML"
                    file_path = os.path.join(folder_mzml, file_name)

                    if not os.path.exists(file_path):
                        new_name_file = f"{sample_name}_0V.mzML"
                        new_file_path = os.path.join(folder_mzml, new_name_file)
                        if os.path.exists(new_file_path):
                            st.warning(f"⚠️ Ficheiro _03 não encontrado! A usar a réplica _02 ({new_name_file}).")
                            file_name = new_name_file
                            file_path = new_file_path

                    if not os.path.exists(file_path):
                        st.error(f"❌ AVISO: Nenhuma réplica encontrada para: {sample_name}")
                        continue

                    st.write(f"📂 A abrir: **{file_name}** -> Analisando {len(sample_info)} targets...")

                    ms1_scans = []
                    with mzml.read(file_path) as reader:
                        for spectra in reader:
                            if spectra["ms level"] == 1:
                                ms1_scans.append({
                                    'rt': spectra['scanList']['scan'][0]['scan start time'],
                                    'mz_array': spectra['m/z array'],
                                    'int_array': spectra['intensity array']
                                })

                    for idx, row in sample_info.iterrows():
                        target_id = row['Target']     
                        target_mz = float(row['mz'])
                        target_rt = float(row['rt'])
                        
                        minor_diference_rt = float('inf') 
                        better_spectra = None
                        
                        for scan in ms1_scans:
                            rt_diference = abs(scan['rt'] - target_rt) 
                            if rt_diference <= tolerance_rt_minutes and rt_diference < minor_diference_rt:
                                minor_diference_rt = rt_diference
                                better_spectra = scan
                        
                        # Criar Report Base
                        report = {
                            "Target": target_id,
                            "Sample": file_name,
                            "Target_m/z": target_mz,
                            "Target_RT": target_rt,
                            "Status": "Not Found",
                            "Main_Intensity": None,
                            "Detected_m/z": None,
                            "Error_ppm": None,
                            "Adducts_Found": "",
                            "Halogen_Diagnosis": "None",
                            "Manual_Isotopic_Envelope": "",
                            "Isotope_Details": "",
                            col_13C: "Not detected",
                            "Estimated_Carbons": None,
                            col_15N: "Not detected",
                            col_18O: "Not detected",
                            col_34S: "Not detected"
                        }

                        if better_spectra:
                            mz_array = better_spectra['mz_array']
                            int_array = better_spectra['int_array']
                            
                            tolerance_mz_da = (target_mz * tolerance_ppm) / 1_000_000
                            valid_indices = np.where((mz_array >= target_mz - tolerance_mz_da) & (mz_array <= target_mz + tolerance_mz_da))[0]
                            
                            if len(valid_indices) > 0:
                                target_index = valid_indices[np.argmax(int_array[valid_indices])]
                                detected_mz = mz_array[target_index]
                                detected_int = int_array[target_index]
                                error_ppm = (abs(detected_mz - target_mz) / target_mz) * 1_000_000
                                
                                report.update({
                                    "Status": "Found",
                                    "Detected_m/z": detected_mz,
                                    "Error_ppm": error_ppm,
                                    "Main_Intensity": detected_int
                                })

                                # --- CARBON-13 ---
                                mz_13c_theoretical = detected_mz + diff_13C
                                tol_13c = (mz_13c_theoretical * tolerance_ppm) / 1_000_000
                                idx_13c = np.where((mz_array >= mz_13c_theoretical - tol_13c) & (mz_array <= mz_13c_theoretical + tol_13c))[0]
                                if len(idx_13c) > 0:
                                    idx_max_13c = idx_13c[np.argmax(int_array[idx_13c])]
                                    abundance_13c = (int_array[idx_max_13c] / detected_int) * 100
                                    diff_13c_real = mz_array[idx_max_13c] - detected_mz
                                    report[col_13C] = f"+{diff_13c_real:.5f} Da ({abundance_13c:.1f}%)"
                                    report["Estimated_Carbons"] = round(abundance_13c / 1.1)

                                # --- HETEROATOMS (N, O, S) ---
                                extra_isotopes = {col_15N: diff_15N, col_18O: diff_18O, col_34S: diff_34S}
                                for iso_name, mass_diff in extra_isotopes.items():
                                    mz_iso_theo = detected_mz + mass_diff
                                    tol_iso = (mz_iso_theo * tolerance_ppm) / 1_000_000
                                    idx_iso = np.where((mz_array >= mz_iso_theo - tol_iso) & (mz_array <= mz_iso_theo + tol_iso))[0]
                                    if len(idx_iso) > 0:
                                        idx_max_iso = idx_iso[np.argmax(int_array[idx_iso])]
                                        abundance_iso = (int_array[idx_max_iso] / detected_int) * 100
                                        diff_real = mz_array[idx_max_iso] - detected_mz
                                        report[iso_name] = f"+{diff_real:.5f} Da ({abundance_iso:.1f}%)"

                                # --- DETAILED ISOTOPE LOGGING ---
                                found_peaks_str = []
                                target_distances = [diff_34S, diff_37Cl, diff_81Br, diff_18O, diff_Cl2, diff_ClBr, diff_Br2, diff_Cl3]
                                registered_mzs = []
                                for dist in target_distances:
                                    mz_iso_target = detected_mz + dist
                                    tol_da = (mz_iso_target * tolerance_ppm) / 1_000_000
                                    idx_found = np.where((mz_array >= mz_iso_target - tol_da) & (mz_array <= mz_iso_target + tol_da))[0]
                                    if len(idx_found) > 0:
                                        idx_max = idx_found[np.argmax(int_array[idx_found])]
                                        mz_found = mz_array[idx_max]
                                        if mz_found not in registered_mzs:
                                            registered_mzs.append(mz_found)
                                            abundance = (int_array[idx_max] / detected_int) * 100
                                            diff_real = mz_found - detected_mz
                                            found_peaks_str.append(f"+{diff_real:.5f} Da ({abundance:.1f}%)")
                                report["Isotope_Details"] = " | ".join(found_peaks_str) if found_peaks_str else "None"

                                # --- MANUAL ENVELOPE ---
                                envelope_str = []
                                idx_env = np.where((mz_array >= detected_mz + 0.5) & (mz_array <= detected_mz + 6.5))[0]
                                for i in idx_env:
                                    abundance_manual = (int_array[i] / detected_int) * 100
                                    if abundance_manual >= 1.0: 
                                        diff_real_manual = mz_array[i] - detected_mz
                                        envelope_str.append(f"+{diff_real_manual:.4f} ({abundance_manual:.1f}%)")
                                report["Manual_Isotopic_Envelope"] = " | ".join(envelope_str) if envelope_str else "No peaks > 1%"

                                # --- HALOGEN DIAGNOSIS ---
                                patterns = {
                                    "1 Chlorine (Cl)": {"mzs": [diff_37Cl], "abundances": [32.0]},
                                    "2 Chlorines (Cl2)": {"mzs": [diff_37Cl, diff_Cl2], "abundances": [65.0, 11.0]},
                                    "3 Chlorines (Cl3)": {"mzs": [diff_37Cl, diff_Cl2, diff_Cl3], "abundances": [98.0, 32.0, 3.0]},
                                    "1 Bromine (Br)": {"mzs": [diff_81Br], "abundances": [98.0]},
                                    "2 Bromines (Br2)": {"mzs": [diff_81Br, diff_Br2], "abundances": [195.0, 95.0]},
                                    "1 Cl + 1 Br": {"mzs": [diff_81Br, diff_ClBr], "abundances": [130.0, 31.0]}
                                }
                                best_match = "None"
                                minor_error = float('inf')
                                for pattern_name, data in patterns.items():
                                    measured_abundances = []
                                    is_valid = True
                                    for diff_mz in data["mzs"]:
                                        mz_target_iso = detected_mz + diff_mz
                                        tol_da = (mz_target_iso * tolerance_ppm) / 1_000_000
                                        idx_enc = np.where((mz_array >= mz_target_iso - tol_da) & (mz_array <= mz_target_iso + tol_da))[0]
                                        
                                        if len(idx_enc) > 0:
                                            idx_max = idx_enc[np.argmax(int_array[idx_enc])]
                                            abundance = (int_array[idx_max] / detected_int) * 100
                                            measured_abundances.append(abundance)
                                        else:
                                            is_valid = False 
                                            break
                                    
                                    if is_valid:
                                        total_error = sum(abs(med - theo) for med, theo in zip(measured_abundances, data["abundances"]))
                                        if total_error < minor_error and total_error <= 20.0:
                                            minor_error = total_error
                                            best_match = pattern_name
                                report["Halogen_Diagnosis"] = best_match

                                # --- ADDUCTS SEARCH ---
                                neutral_mass = detected_mz - 1.007276 
                                adducts_str = []
                                for _, r_ad in df_adducts.iterrows():
                                    if r_ad['Ion name'] == "[M+H]+": continue
                                    theoretical_mz = (neutral_mass * r_ad['Mult']) + r_ad['Mass']
                                    tol_adduct = (theoretical_mz * tolerance_ppm) / 1_000_000
                                    idx_ad = np.where((mz_array >= theoretical_mz - tol_adduct) & (mz_array <= theoretical_mz + tol_adduct))[0]
                                    if len(idx_ad) > 0: 
                                        adducts_str.append(r_ad['Ion name'])
                                report["Adducts_Found"] = ", ".join(adducts_str)

                        final_results.append(report)

                    # Limpar memória
                    del ms1_scans
                    gc.collect()
                    
                    # Atualizar barra de progresso
                    progress_bar.progress((idx_sample + 1) / total_samples)
                
                status.update(label="Análise concluída!", state="complete", expanded=False)

            except Exception as e:
                st.error(f"Ocorreu um erro durante a execução: {e}")
                status.update(label="Erro na análise", state="error")

        # ==========================================
        # 3. Disponibilizar Excel para Download
        # ==========================================
        if 'final_results' in locals() and final_results:
            st.success("✅ Análise finalizada com sucesso! Podes descarregar o ficheiro abaixo.")
            
            df_final = pd.DataFrame(final_results)
            
            # Criar o Excel em memória (RAM) em vez de guardar no disco
            output = BytesIO()
            with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
                df_final.to_excel(writer, index=False, sheet_name="Resultados")
            output.seek(0)
            
            st.download_button(
                label="📥 Descarregar Resultados (Excel)",
                data=output,
                file_name="resultados_metabolomica.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )