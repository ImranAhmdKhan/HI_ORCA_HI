import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import re
import math # For exp in Curtin-Hammett

# Try to import matplotlib
try:
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator # For better integer ticks on x-axis if needed
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

# --- Conversion Factors ---
HARTREE_TO_KCAL_PER_MOL = 627.50960803
HARTREE_TO_KJ_PER_MOL = 2625.49963948
HARTREE_TO_EV = 27.21139664
HARTREE_TO_J = 4.3597447222071e-18
AVOGADRO_CONSTANT = 6.02214076e23
J_PER_MOL_K_FROM_HARTREE_PER_K = HARTREE_TO_J * AVOGADRO_CONSTANT
CAL_PER_MOL_K_FROM_J_PER_MOL_K = 1 / 4.184
R_KCAL_MOL_K = 1.987204259e-3 # kcal/mol·K


def parse_orca_thermo_block(filepath):
    # (This function remains identical to the previous version with XYZ fallback)
    thermo_data = {
        "temperature_k": None, "electronic_energy_eh": None, "zero_point_energy_eh": None,
        "thermal_vibrational_correction_eh": None, "thermal_rotational_correction_eh": None,
        "thermal_translational_correction_eh": None, "total_inner_energy_u_eh": None,
        "thermal_enthalpy_correction_eh": None, "total_enthalpy_h_eh": None,
        "ts_electronic_eh": None, "ts_vibrational_eh": None, "ts_rotational_eh": None,
        "ts_translational_eh": None, "final_ts_term_eh": None,
        "final_gibbs_free_energy_g_eh": None, "g_minus_eel_eh": None,
        "total_thermal_correction_to_e_eh": None, "total_correction_to_e_eh": None,
        "source_info": f"Not Parsed ({os.path.basename(filepath)})", "filepath": filepath
    }
    patterns = {
        "temperature_k": re.compile(r"THERMOCHEMISTRY AT\s+([\d\.]+)\s*K"),
        "electronic_energy_eh": re.compile(r"Electronic energy\s+\.+\s+(-?[\d\.]+)\s*Eh"),
        "zero_point_energy_eh": re.compile(r"Zero point energy\s+\.+\s+(-?[\d\.]+)\s*Eh"),
        "thermal_vibrational_correction_eh": re.compile(r"Thermal vibrational correction\s+\.+\s+(-?[\d\.]+)\s*Eh"),
        "thermal_rotational_correction_eh": re.compile(r"Thermal rotational correction\s+\.+\s+(-?[\d\.]+)\s*Eh"),
        "thermal_translational_correction_eh": re.compile(r"Thermal translational correction\s+\.+\s+(-?[\d\.]+)\s*Eh"),
        "total_inner_energy_u_eh": re.compile(r"Total thermal energy\s+\.+\s+(-?[\d\.]+)\s*Eh"), # ORCA < 5 might use "Inner energy"
        "total_thermal_correction_to_e_eh": re.compile(r"Total thermal correction\s+\.+\s+(-?[\d\.]+)\s*Eh"),
        "total_correction_to_e_eh": re.compile(r"Total correction\s+\.+\s+(-?[\d\.]+)\s*Eh"), # ZPE + thermal
        "thermal_enthalpy_correction_eh": re.compile(r"Thermal Enthalpy correction\s+\.+\s+(-?[\d\.]+)\s*Eh"), # This is kBT for ideal gas
        "total_enthalpy_h_eh": re.compile(r"Total Enthalpy\s+\.+\s+(-?[\d\.]+)\s*Eh"), # H = U + kBT
        "ts_electronic_eh": re.compile(r"Electronic entropy\s+\.+\s+(-?[\d\.]+)\s*Eh"),
        "ts_vibrational_eh": re.compile(r"Vibrational entropy\s+\.+\s+(-?[\d\.]+)\s*Eh"),
        "ts_rotational_eh": re.compile(r"Rotational entropy\s+\.+\s+(-?[\d\.]+)\s*Eh"),
        "ts_translational_eh": re.compile(r"Translational entropy\s+\.+\s+(-?[\d\.]+)\s*Eh"),
        "final_ts_term_eh": re.compile(r"Final entropy term\s+\.+\s+(-?[\d\.]+)\s*Eh"), # This is T*S_total
        "final_gibbs_free_energy_g_eh": re.compile(r"Final Gibbs free energy\s+\.+\s+(-?[\d\.]+)\s*Eh"), # G = H - T*S_total
        "g_minus_eel_eh": re.compile(r"G-E\(el\)\s+\.+\s+(-?[\d\.]+)\s*Eh") # G - E_el
    }
    found_thermo_block_start = False; parsed_from_thermo_block = False
    lines = []
    try:
        try:
            with open(filepath, 'r', encoding='utf-8') as f: lines = f.readlines()
        except UnicodeDecodeError:
            print(f"Warning: UTF-8 decoding failed for {os.path.basename(filepath)}. Trying latin-1.")
            with open(filepath, 'r', encoding='latin-1') as f: lines = f.readlines()
        
        for i, line in enumerate(lines):
            if not found_thermo_block_start:
                temp_match = patterns["temperature_k"].search(line)
                if temp_match:
                    thermo_data["temperature_k"] = float(temp_match.group(1)); found_thermo_block_start = True
                    parsed_from_thermo_block = True; thermo_data["source_info"] = f"Full Thermo ({os.path.basename(filepath)})"; continue
            if found_thermo_block_start:
                # Attempt to parse "Inner energy" for older ORCA versions if "Total thermal energy" isn't found yet
                if thermo_data.get("total_inner_energy_u_eh") is None:
                    inner_e_match = re.search(r"Inner energy\s+\.+\s+(-?[\d\.]+)\s*Eh", line)
                    if inner_e_match:
                        thermo_data["total_inner_energy_u_eh"] = float(inner_e_match.group(1))

                for key, pattern in patterns.items():
                    if key == "temperature_k": continue
                    match = pattern.search(line)
                    if match:
                        try: thermo_data[key] = float(match.group(1))
                        except ValueError: pass 
                if "G-E(el)" in line and thermo_data["g_minus_eel_eh"] is not None: # Check if we are at the end of thermo block
                    look_ahead_lines = lines[i+1 : i+6]; 
                    if any("SUGGESTED CITATIONS" in nl for nl in look_ahead_lines) or \
                       any("Timings for individual modules" in nl for nl in look_ahead_lines) or \
                       any("Total run time" in nl for nl in look_ahead_lines): # Added more end markers
                        break 
        
        if not parsed_from_thermo_block or thermo_data.get("total_inner_energy_u_eh") is None : # If full block not parsed or U is missing
            xyz_filepath = os.path.splitext(filepath)[0] + ".xyz"; xyz_energy_parsed = False
            if os.path.exists(xyz_filepath):
                xyz_lines = []
                try:
                    try:
                        with open(xyz_filepath, 'r', encoding='utf-8') as xyz_f: xyz_lines = xyz_f.readlines()
                    except UnicodeDecodeError:
                        print(f"Warning: UTF-8 decoding failed for XYZ {os.path.basename(xyz_filepath)}. Trying latin-1.")
                        with open(xyz_filepath, 'r', encoding='latin-1') as xyz_f: xyz_lines = xyz_f.readlines()
                    
                    if len(xyz_lines) >= 2:
                        comment_line = xyz_lines[1].strip()
                        # More robust regex for energy in XYZ comment:
                        energy_match = re.search(r"(?:Energy|E|energy|E\(SCF\)|SCF)\s*=?\s*(-?\d+\.?\d*(?:[eE][-+]?\d+)?)", comment_line, re.IGNORECASE)
                        if energy_match:
                            xyz_energy = float(energy_match.group(1))
                            if not parsed_from_thermo_block: # Only use XYZ if no thermo block was found
                                thermo_data["electronic_energy_eh"] = xyz_energy
                                thermo_data["source_info"] = f"XYZ Comment ({os.path.basename(xyz_filepath)})"; xyz_energy_parsed = True
                            # If thermo block was parsed but some values are missing, we don't overwrite E_el from XYZ
                            # but we note that XYZ energy was found.
                            elif thermo_data.get("electronic_energy_eh") is None: # If E_el is missing from incomplete thermo
                                thermo_data["electronic_energy_eh"] = xyz_energy
                                # Keep source_info as "Full Thermo" but it's incomplete
                                print(f"Note: Used E_el from XYZ for incomplete thermo block of {os.path.basename(filepath)}")
                                xyz_energy_parsed = True


                except Exception as e_xyz: print(f"Warning: Could not parse XYZ {xyz_filepath}: {e_xyz}")
            
            if not parsed_from_thermo_block and not xyz_energy_parsed: 
                messagebox.showinfo("Parsing Info", f"Neither full thermochemistry nor XYZ comment energy found for '{os.path.basename(filepath)}'."); return None
            elif parsed_from_thermo_block and thermo_data.get("total_inner_energy_u_eh") is None and thermo_data.get("electronic_energy_eh") is None:
                 messagebox.showwarning("Parsing Warning", f"Incomplete thermochemistry block in '{os.path.basename(filepath)}' and no E_el. Some values may be missing or N/A.")
            elif parsed_from_thermo_block and thermo_data.get("total_inner_energy_u_eh") is None:
                 print(f"Note: Incomplete thermochemistry block in '{os.path.basename(filepath)}'. Some values may be missing.")


    except FileNotFoundError: messagebox.showerror("Error", f"File not found: {filepath}"); return None
    except Exception as e: messagebox.showerror("Error", f"Error parsing '{os.path.basename(filepath)}':\n{e}"); import traceback; traceback.print_exc(); return None
    return thermo_data


def populate_treeview_widget(treeview_widget, data, app_instance):
    # (This function remains identical to the previous version for displaying single file thermo data)
    for item in treeview_widget.get_children(): treeview_widget.delete(item)
    if data is None: treeview_widget.insert('', tk.END, values=("Error", "No data for selected file.", "", "")); return

    def format_num(value, precision=8, is_scientific=False):
        if value is None: return "N/A"
        return f"{value:.{precision}e}" if is_scientific else f"{value:.{precision}f}"

    source_info_str = data.get("source_info", ""); is_xyz_only_data = "XYZ Comment" in source_info_str

    if app_instance.show_general_var.get():
        category_gen = "General Information"; temp_val = data.get('temperature_k')
        if temp_val is not None and not is_xyz_only_data: treeview_widget.insert('', tk.END, values=(category_gen, "Temperature", format_num(temp_val, 2), "K"))
        treeview_widget.insert('', tk.END, values=(category_gen, "Data Source", source_info_str, ""))
        treeview_widget.insert('', tk.END, values=("---", "---", "---", "---"), tags=('separator',))

    def insert_energy_param_tree(cat, param_base, val_eh, is_sum=False):
        param_display_name_base = f"{param_base} [SUM]" if is_sum else param_base
        if val_eh is None: treeview_widget.insert('', tk.END, values=(cat, param_display_name_base, "N/A", "Eh"))
        else:
            units_to_display_data = []
            if app_instance.show_eh_var.get(): units_to_display_data.append({'val': val_eh, 'unit': "Eh", 'prec': 8, 'is_sci': False})
            if app_instance.show_kcal_var.get(): units_to_display_data.append({'val': val_eh * HARTREE_TO_KCAL_PER_MOL, 'unit': "kcal/mol", 'prec': 2, 'is_sci': False})
            if app_instance.show_kj_var.get(): units_to_display_data.append({'val': val_eh * HARTREE_TO_KJ_PER_MOL, 'unit': "kJ/mol", 'prec': 2, 'is_sci': False})
            if app_instance.show_ev_var.get(): units_to_display_data.append({'val': val_eh * HARTREE_TO_EV, 'unit': "eV", 'prec': 4, 'is_sci': False})
            if not units_to_display_data: treeview_widget.insert('', tk.END, values=(cat, param_display_name_base, " (No units selected)", ""))
            else:
                for i, unit_data in enumerate(units_to_display_data):
                    display_name_row = param_display_name_base if i == 0 and len(units_to_display_data) == 1 else \
                                       (param_display_name_base if i == 0 else f"  └─ ({unit_data['unit']})")
                    if len(units_to_display_data) == 1: display_name_row = param_display_name_base
                    treeview_widget.insert('', tk.END, values=(cat if i == 0 else "", display_name_row,
                                           format_num(unit_data['val'], unit_data['prec'], unit_data['is_sci']), unit_data['unit']))
        if is_sum : treeview_widget.insert('', tk.END, values=("---", "---", "---", "---"), tags=('separator',))

    if app_instance.show_inner_energy_var.get():
        category = "Inner Energy (U)"; eel_param_name = "Electronic energy (E_el)"
        insert_energy_param_tree(category, eel_param_name, data.get("electronic_energy_eh"))
        if not is_xyz_only_data:
            insert_energy_param_tree(category, "Zero point energy (E_ZPE)", data.get("zero_point_energy_eh"))
            insert_energy_param_tree(category, "Thermal vibrational correction (E_vib)", data.get("thermal_vibrational_correction_eh"))
            insert_energy_param_tree(category, "Thermal rotational correction (E_rot)", data.get("thermal_rotational_correction_eh"))
            insert_energy_param_tree(category, "Thermal translational correction (E_trans)", data.get("thermal_translational_correction_eh"))
            insert_energy_param_tree(category, "Total Inner Energy (U)", data.get("total_inner_energy_u_eh"), is_sum=True)
        elif data.get("electronic_energy_eh") is not None :
             insert_energy_param_tree(category, "Total Inner Energy (U) (approx. as E_el)", data.get("electronic_energy_eh"), is_sum=True)
        else: treeview_widget.insert('', tk.END, values=(category, "Total Inner Energy (U)", "N/A", "")); treeview_widget.insert('', tk.END, values=("---", "---", "---", "---"), tags=('separator',))

    sections_to_display = [
        (app_instance.show_corrections_var, "Overall Energy Corrections (to E_el)", [
            ("Total thermal correction to E_el", "total_thermal_correction_to_e_eh", False),
            ("Total correction (ZPE + thermal)", "total_correction_to_e_eh", True)]),
        (app_instance.show_enthalpy_var, "Enthalpy (H)", [
            ("Thermal Enthalpy corr. (kB*T term)", "thermal_enthalpy_correction_eh", False),
            ("Total Enthalpy (H)", "total_enthalpy_h_eh", True)]),
        (app_instance.show_ts_terms_var, "Entropy (T*S terms)", [
            ("T*S (electronic)", "ts_electronic_eh", False), ("T*S (vibrational)", "ts_vibrational_eh", False),
            ("T*S (rotational)", "ts_rotational_eh", False), ("T*S (translational)", "ts_translational_eh", False),
            ("Final T*S term (Total T*S)", "final_ts_term_eh", True)]),
        (app_instance.show_gibbs_var, "Gibbs Free Energy (G)", [
            ("Total -T*S term (entropy correction)", -data.get("final_ts_term_eh") if data.get("final_ts_term_eh") is not None else None, False),
            ("Final Gibbs Free Energy (G)", "final_gibbs_free_energy_g_eh", True),
            ("G - E_el", "g_minus_eel_eh", True)])]

    if not is_xyz_only_data:
        for show_var, category, params in sections_to_display:
            if show_var.get():
                for param_name, data_key, is_sum_flag in params:
                    value = data.get(data_key) if isinstance(data_key, str) else data_key
                    insert_energy_param_tree(category, param_name, value, is_sum=is_sum_flag)
    elif is_xyz_only_data: # Only E_el is available
        eel_val = data.get("electronic_energy_eh")
        if app_instance.show_enthalpy_var.get(): insert_energy_param_tree("Enthalpy (H)", "Total Enthalpy (H) (approx. as E_el)", eel_val, is_sum=True)
        if app_instance.show_gibbs_var.get(): insert_energy_param_tree("Gibbs Free Energy (G)", "Final Gibbs Free Energy (G) (approx. as E_el)", eel_val, is_sum=True)

    if app_instance.show_s_total_var.get() and not is_xyz_only_data:
        category = "Total Entropy (S_total)"
        if data.get("final_ts_term_eh") is not None and data.get("temperature_k") is not None and data["temperature_k"] > 0:
            ts_total_eh = data["final_ts_term_eh"]; temp_k = data["temperature_k"]
            s_total_au = ts_total_eh / temp_k; s_j_mol_k = s_total_au * J_PER_MOL_K_FROM_HARTREE_PER_K; s_cal_mol_k = s_j_mol_k * CAL_PER_MOL_K_FROM_J_PER_MOL_K
            treeview_widget.insert('', tk.END, values=(category, "S_total", format_num(s_total_au, 8, True), "Hartree/K"))
            treeview_widget.insert('', tk.END, values=("", f"  └─ (J/mol·K)", format_num(s_j_mol_k, 2), "J/mol·K"))
            treeview_widget.insert('', tk.END, values=("", f"  └─ (cal/mol·K)", format_num(s_cal_mol_k, 2), "cal/mol·K"))
        else: treeview_widget.insert('', tk.END, values=(category, "S_total", "N/A", ""))
        treeview_widget.insert('', tk.END, values=("---", "---", "---", "---"), tags=('separator',))
    app_instance.apply_treeview_row_tags(treeview_widget)


def parse_multiframe_xyz(filepath):
    # (This function remains identical to the previous version)
    frames_data = []
    lines = []
    try:
        try:
            with open(filepath, 'r', encoding='utf-8') as f: lines = f.readlines()
        except UnicodeDecodeError:
            print(f"Warning: UTF-8 decoding failed for multi-frame XYZ {os.path.basename(filepath)}. Trying latin-1.")
            with open(filepath, 'r', encoding='latin-1') as f: lines = f.readlines()

        idx = 0; frame_index = 0
        while idx < len(lines):
            try:
                num_atoms_str = lines[idx].strip(); 
                if not num_atoms_str: idx += 1; continue # Skip empty lines
                num_atoms = int(num_atoms_str); idx += 1
                if idx + num_atoms + 1 > len(lines) +1 : break # Check if enough lines for comment + coords
                comment_line = lines[idx].strip(); idx += 1
                coordinates = [lines[idx+k].strip() for k in range(num_atoms)]; idx += num_atoms
                energy_eh = None
                # More robust regex for energy in XYZ comment:
                energy_match = re.search(r"(?:Energy|E|energy|E\(SCF\)|SCF)\s*=?\s*(-?\d+\.?\d*(?:[eE][-+]?\d+)?)", comment_line, re.IGNORECASE)
                if energy_match:
                    try: energy_eh = float(energy_match.group(1))
                    except ValueError: pass
                frames_data.append({'frame_index': frame_index, 'num_atoms': num_atoms, 'comment_line': comment_line, 
                                    'energy_eh': energy_eh, 'coordinates': coordinates, 'filepath': filepath})
                frame_index += 1
            except ValueError: print(f"Warning: Skipping invalid frame in {os.path.basename(filepath)} near line {idx+1}."); break 
            except IndexError: print(f"Warning: File {os.path.basename(filepath)} ended unexpectedly during frame parsing."); break
    except FileNotFoundError: messagebox.showerror("Error", f"File not found: {filepath}"); return None
    except Exception as e: messagebox.showerror("Error", f"Error reading {os.path.basename(filepath)}: {e}"); return None
    if not frames_data: messagebox.showinfo("Parsing Info", f"No XYZ frames in {os.path.basename(filepath)}."); return []
    return frames_data


class ThermoApp:
    def __init__(self, master):
        self.master = master
        master.title("Hi_ORCA Thermo & TS Tools - File Reorder_Hi(CHEMXELENT)")
        master.geometry("1250x900")

        # Common data storage
        self.loaded_file_paths = []
        self.loaded_files_data = {}

        # Tab 1 specific
        self.current_selected_filepath = None
        self.ch_reactant_var = tk.StringVar()
        self.ch_product_var = tk.StringVar()
        self.plot_reference_var = tk.StringVar()
        self.plot_energy_type_var = tk.StringVar()
        # For consolidated table
        self.consolidated_table_window = None 
        self.consolidated_param_var = tk.StringVar()
        self.consolidated_unit_var = tk.StringVar()
        # For axis break (Tab 1 plot)
        self.enable_axis_break_var = tk.BooleanVar(value=False)
        self.axis_break_from_var = tk.StringVar()
        self.axis_break_to_var = tk.StringVar()
        # For Curtin-Hammett plot axis break
        self.enable_ch_axis_break_var = tk.BooleanVar(value=False)
        self.ch_axis_break_from_var = tk.StringVar()
        self.ch_axis_break_to_var = tk.StringVar()
        # For Step Delta Calculation Window
        self.step_delta_window = None
        self.step_reactant_var = tk.StringVar()
        self.step_product_var = tk.StringVar()
        self.step_delta_results_data = [] # To store [(r_name, p_name, dG, dH, dS, TdS, Tr, Tp), ...]
        # Plot style options for Tab 1 general plot
        self.plot_line_style_var = tk.StringVar(value='-') # solid line
        self.plot_line_width_var = tk.StringVar(value='2.0')
        self.plot_line_color_var = tk.StringVar(value='dodgerblue')
        # Plot style options for CH plot
        self.ch_plot_line_style_var = tk.StringVar(value='-')
        self.ch_plot_line_width_var = tk.StringVar(value='2.0')
        self.ch_plot_line_color_var = tk.StringVar(value='purple')


        # Tab 2 specific
        self.allxyz_file_paths = []
        self.allxyz_files_data = {}
        self.current_selected_allxyz_filepath = None
        self.highest_energy_frame_data = None
        self.enable_allxyz_axis_break_var = tk.BooleanVar(value=False)
        self.allxyz_axis_break_from_var = tk.StringVar()
        self.allxyz_axis_break_to_var = tk.StringVar()
        # Plot style options for Tab 2 allxyz plot
        self.allxyz_plot_line_style_var = tk.StringVar(value='-')
        self.allxyz_plot_line_width_var = tk.StringVar(value='1.5')
        self.allxyz_plot_line_color_var = tk.StringVar(value='green')


        # Common display options for Tab 1
        self.show_eh_var = tk.BooleanVar(value=True); self.show_kcal_var = tk.BooleanVar(value=True)
        self.show_kj_var = tk.BooleanVar(value=True); self.show_ev_var = tk.BooleanVar(value=True)
        self.show_general_var = tk.BooleanVar(value=True); self.show_inner_energy_var = tk.BooleanVar(value=True)
        self.show_corrections_var = tk.BooleanVar(value=True); self.show_enthalpy_var = tk.BooleanVar(value=True)
        self.show_ts_terms_var = tk.BooleanVar(value=True); self.show_s_total_var = tk.BooleanVar(value=True)
        self.show_gibbs_var = tk.BooleanVar(value=True)

        self.row_tags = ('evenrow', 'oddrow')
        self.line_styles = ['-', '--', '-.', ':'] # solid, dashed, dashdot, dotted
        self.plot_colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', 'purple', 'orange', 'brown', 'dodgerblue']


        style = ttk.Style(master)
        style.configure("Treeview.Heading", font=('Calibri', 10, 'bold'))
        style.map('Treeview', background=[('selected', '#BFBFBF')])

        self.notebook = ttk.Notebook(master)
        self.notebook.pack(expand=True, fill='both', padx=5, pady=5)

        self.tab1 = ttk.Frame(self.notebook, padding=5)
        self.notebook.add(self.tab1, text='Thermochemistry & Reaction Analysis')
        self.setup_tab1_ui()

        self.tab2 = ttk.Frame(self.notebook, padding=5)
        self.notebook.add(self.tab2, text='TS Estimation (.allxyz)')
        self.setup_tab2_ui()

    # --- Helper methods for energy retrieval ---
    def _get_relevant_energy(self, data_dict, energy_type_priority=['G', 'H', 'E_el+ZPE', 'E_el']):
        if not data_dict: return None, "N/A", None
        temp_k = data_dict.get("temperature_k", 298.15) # Default temp if not found
        if temp_k is None: temp_k = 298.15 # Ensure temp_k is not None

        for etype in energy_type_priority:
            if etype == 'G' and data_dict.get("final_gibbs_free_energy_g_eh") is not None:
                return data_dict["final_gibbs_free_energy_g_eh"], "G", temp_k
            if etype == 'H' and data_dict.get("total_enthalpy_h_eh") is not None:
                return data_dict["total_enthalpy_h_eh"], "H", temp_k
            if etype == 'E_el+ZPE' and data_dict.get("electronic_energy_eh") is not None and data_dict.get("zero_point_energy_eh") is not None:
                return data_dict["electronic_energy_eh"] + data_dict["zero_point_energy_eh"], "E+ZPE", temp_k
            if etype == 'E_el' and data_dict.get("electronic_energy_eh") is not None:
                return data_dict["electronic_energy_eh"], "E_el", temp_k
        return None, "N/A", temp_k

    def _get_s_and_ts(self, data_dict):
        if data_dict and data_dict.get("final_ts_term_eh") is not None and data_dict.get("temperature_k") is not None:
            ts_eh = data_dict["final_ts_term_eh"]; temp = data_dict["temperature_k"]
            if temp is not None and temp > 0: return ts_eh / temp, ts_eh, temp
        default_temp = data_dict.get("temperature_k", 298.15)
        if default_temp is None: default_temp = 298.15
        return None, None, default_temp


    def _get_energy_key_for_plot(self, requested_type_str, data_dict):
        if not data_dict: return None
        if requested_type_str == "Gibbs Free Energy (G)": return "final_gibbs_free_energy_g_eh"
        if requested_type_str == "Enthalpy (H)": return "total_enthalpy_h_eh"
        if requested_type_str == "E_el + ZPVE":
            # This will be handled by directly calculating E_el + ZPVE if available
            if data_dict.get("electronic_energy_eh") is not None and data_dict.get("zero_point_energy_eh") is not None:
                return "_e_plus_zpe_calculated" # Special key, not directly in data_dict
            return "electronic_energy_eh" # Fallback if ZPVE is missing
        if requested_type_str == "Electronic Energy (E_el)": return "electronic_energy_eh"
        return "electronic_energy_eh" # Default fallback


    def setup_tab1_ui(self):
        master = self.tab1
        self.main_paned_window = tk.PanedWindow(master, orient=tk.HORIZONTAL, sashrelief=tk.RAISED)
        self.main_paned_window.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        self.left_pane = ttk.Frame(self.main_paned_window, padding=5)
        self.main_paned_window.add(self.left_pane, width=350, minsize=250)

        tk.Button(self.left_pane, text="Load ORCA File(s) (.out/.log/.xyz)", command=self.load_files).pack(pady=(0,5), fill=tk.X)

        listbox_control_frame = ttk.Frame(self.left_pane)
        listbox_control_frame.pack(pady=(0,5), fill=tk.BOTH, expand=True)

        listbox_frame = ttk.Frame(listbox_control_frame); listbox_frame.pack(pady=0, fill=tk.BOTH, expand=True)
        self.file_listbox = tk.Listbox(listbox_frame, selectmode=tk.SINGLE, exportselection=False, height=8)
        self.file_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        listbox_scrollbar = ttk.Scrollbar(listbox_frame, orient="vertical", command=self.file_listbox.yview)
        listbox_scrollbar.pack(side=tk.RIGHT, fill=tk.Y); self.file_listbox.config(yscrollcommand=listbox_scrollbar.set)
        self.file_listbox.bind('<<ListboxSelect>>', self.on_file_select_from_listbox)

        reorder_buttons_frame = ttk.Frame(listbox_control_frame)
        reorder_buttons_frame.pack(fill=tk.X, pady=(5,0))
        ttk.Button(reorder_buttons_frame, text="Move Up ↑", command=lambda: self.move_file_in_list("up")).pack(side=tk.LEFT, expand=True, padx=2)
        ttk.Button(reorder_buttons_frame, text="Move Down ↓", command=lambda: self.move_file_in_list("down")).pack(side=tk.LEFT, expand=True, padx=2)

        options_outer_frame = ttk.Frame(self.left_pane); options_outer_frame.pack(fill=tk.X, pady=5)
        self.unit_frame = ttk.LabelFrame(options_outer_frame, text="Energy Units (Table Display)", padding=(10, 5))
        self.unit_frame.pack(pady=5, fill=tk.X)
        ttk.Checkbutton(self.unit_frame, text="Eh", variable=self.show_eh_var, command=self.on_selection_change).grid(row=0, column=0, sticky='w', padx=2)
        ttk.Checkbutton(self.unit_frame, text="kcal/mol", variable=self.show_kcal_var, command=self.on_selection_change).grid(row=0, column=1, sticky='w', padx=2)
        ttk.Checkbutton(self.unit_frame, text="kJ/mol", variable=self.show_kj_var, command=self.on_selection_change).grid(row=1, column=0, sticky='w', padx=2)
        ttk.Checkbutton(self.unit_frame, text="eV", variable=self.show_ev_var, command=self.on_selection_change).grid(row=1, column=1, sticky='w', padx=2)

        self.section_frame = ttk.LabelFrame(options_outer_frame, text="Thermo Sections (Table Display)", padding=(10, 5))
        self.section_frame.pack(pady=5, fill=tk.X)
        sec_cb_frame1 = ttk.Frame(self.section_frame); sec_cb_frame1.pack(side=tk.LEFT, anchor='nw', padx=2)
        sec_cb_frame2 = ttk.Frame(self.section_frame); sec_cb_frame2.pack(side=tk.LEFT, anchor='nw', padx=2)
        ttk.Checkbutton(sec_cb_frame1, text="General", variable=self.show_general_var, command=self.on_selection_change).pack(anchor='w')
        ttk.Checkbutton(sec_cb_frame1, text="Inner Energy", variable=self.show_inner_energy_var, command=self.on_selection_change).pack(anchor='w')
        ttk.Checkbutton(sec_cb_frame1, text="E_el Correct.", variable=self.show_corrections_var, command=self.on_selection_change).pack(anchor='w')
        ttk.Checkbutton(sec_cb_frame1, text="Enthalpy", variable=self.show_enthalpy_var, command=self.on_selection_change).pack(anchor='w')
        ttk.Checkbutton(sec_cb_frame2, text="T*S Terms", variable=self.show_ts_terms_var, command=self.on_selection_change).pack(anchor='w')
        ttk.Checkbutton(sec_cb_frame2, text="Total Entropy", variable=self.show_s_total_var, command=self.on_selection_change).pack(anchor='w')
        ttk.Checkbutton(sec_cb_frame2, text="Gibbs Energy", variable=self.show_gibbs_var, command=self.on_selection_change).pack(anchor='w')

        self.right_pane = ttk.Frame(self.main_paned_window, padding=5)
        self.main_paned_window.add(self.right_pane)
        self.filepath_display_frame = ttk.Frame(self.right_pane)
        self.filepath_display_frame.pack(fill=tk.X, pady=(0,5))
        self.filepath_label = tk.Label(self.filepath_display_frame, text="No file selected for display", relief=tk.SUNKEN, width=50, anchor='w')
        self.filepath_label.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.display_button = tk.Button(self.filepath_display_frame, text="Refresh Display", command=self.trigger_refresh_display, state=tk.DISABLED)
        self.display_button.pack(side=tk.LEFT, padx=(10,0))

        action_buttons_frame = ttk.Frame(self.right_pane); action_buttons_frame.pack(fill=tk.X, pady=5)
        self.copy_button = ttk.Button(action_buttons_frame, text="Copy Thermo Table (TSV)", command=self.copy_thermo_table_to_clipboard)
        self.copy_button.pack(side=tk.LEFT, padx=(0,5))

        # Consolidated Table Button
        self.show_all_results_button = ttk.Button(action_buttons_frame, text="Consolidated Table", command=self.open_consolidated_table_window)
        self.show_all_results_button.pack(side=tk.LEFT, padx=(5,0))
        # Step Deltas Button
        self.calculate_step_deltas_button = ttk.Button(action_buttons_frame, text="Step Δ Values", command=self.open_step_deltas_window)
        self.calculate_step_deltas_button.pack(side=tk.LEFT, padx=(5,0))


        plot_options_outer_frame = ttk.LabelFrame(self.right_pane, text="Plotting Options (General Profile)", padding=10)
        plot_options_outer_frame.pack(fill=tk.X, pady=5)

        plot_options_frame = ttk.Frame(plot_options_outer_frame) # For Ref and Energy Type
        plot_options_frame.pack(fill=tk.X)

        ttk.Label(plot_options_frame, text="Plot Ref:").grid(row=0, column=0, padx=(0,2), pady=2, sticky=tk.W)
        self.plot_reference_combo = ttk.Combobox(plot_options_frame, textvariable=self.plot_reference_var, state="readonly", width=20)
        self.plot_reference_combo.grid(row=0, column=1, padx=(0,5), pady=2, sticky=tk.EW)
        
        ttk.Label(plot_options_frame, text="Plot Energy Type:").grid(row=0, column=2, padx=(5,2), pady=2, sticky=tk.W)
        self.plot_energy_type_combo = ttk.Combobox(plot_options_frame, textvariable=self.plot_energy_type_var,
                                                values=["Gibbs Free Energy (G)", "Enthalpy (H)", "E_el + ZPVE", "Electronic Energy (E_el)"],
                                                state="readonly", width=22)
        self.plot_energy_type_combo.grid(row=0, column=3, padx=(0,5), pady=2, sticky=tk.EW)
        self.plot_energy_type_var.set("Gibbs Free Energy (G)")
        self.plot_energy_type_combo.bind("<<ComboboxSelected>>", lambda e: self.check_plot_button_state())
        
        self.plot_energy_button = ttk.Button(plot_options_frame, text="Plot Profile", command=self.plot_energy_profile, state=tk.DISABLED)
        self.plot_energy_button.grid(row=0, column=4, padx=(5,0), pady=2, sticky=tk.E)
        
        plot_options_frame.grid_columnconfigure(1, weight=1)
        plot_options_frame.grid_columnconfigure(3, weight=1)

        # --- Plot Style Options for Tab 1 General Plot ---
        plot_style_frame = ttk.Frame(plot_options_outer_frame)
        plot_style_frame.pack(fill=tk.X, pady=(5,0))
        ttk.Label(plot_style_frame, text="Line Style:").grid(row=0, column=0, padx=(0,2), pady=2, sticky=tk.W)
        self.plot_line_style_combo = ttk.Combobox(plot_style_frame, textvariable=self.plot_line_style_var, values=self.line_styles, state="readonly", width=8)
        self.plot_line_style_combo.grid(row=0, column=1, padx=(0,5), pady=2, sticky=tk.W)
        ttk.Label(plot_style_frame, text="Line Width:").grid(row=0, column=2, padx=(5,2), pady=2, sticky=tk.W)
        self.plot_line_width_entry = ttk.Entry(plot_style_frame, textvariable=self.plot_line_width_var, width=5)
        self.plot_line_width_entry.grid(row=0, column=3, padx=(0,5), pady=2, sticky=tk.W)
        ttk.Label(plot_style_frame, text="Color:").grid(row=0, column=4, padx=(5,2), pady=2, sticky=tk.W)
        self.plot_line_color_combo = ttk.Combobox(plot_style_frame, textvariable=self.plot_line_color_var, values=self.plot_colors, state="readonly", width=10)
        self.plot_line_color_combo.grid(row=0, column=5, padx=(0,5), pady=2, sticky=tk.W)


        # --- Axis Break UI Elements for Tab 1 ---
        axis_break_frame = ttk.Frame(plot_options_outer_frame)
        axis_break_frame.pack(fill=tk.X, pady=(5,0))

        self.enable_axis_break_cb = ttk.Checkbutton(axis_break_frame, text="Enable Y-Axis Break", variable=self.enable_axis_break_var, command=self._toggle_axis_break_entries_tab1)
        self.enable_axis_break_cb.grid(row=0, column=0, padx=(0,5), pady=2, sticky=tk.W)

        ttk.Label(axis_break_frame, text="Break From (kcal/mol):").grid(row=0, column=1, padx=(0,2), pady=2, sticky=tk.W)
        self.axis_break_from_entry = ttk.Entry(axis_break_frame, textvariable=self.axis_break_from_var, width=8, state=tk.DISABLED)
        self.axis_break_from_entry.grid(row=0, column=2, padx=(0,5), pady=2, sticky=tk.W)

        ttk.Label(axis_break_frame, text="Break To (kcal/mol):").grid(row=0, column=3, padx=(0,2), pady=2, sticky=tk.W)
        self.axis_break_to_entry = ttk.Entry(axis_break_frame, textvariable=self.axis_break_to_var, width=8, state=tk.DISABLED)
        self.axis_break_to_entry.grid(row=0, column=4, padx=(0,0), pady=2, sticky=tk.W)
        # --- End Axis Break UI ---


        tree_frame = ttk.Frame(self.right_pane); tree_frame.pack(fill=tk.BOTH, expand=True, pady=(5,0))
        columns = ("category", "parameter", "value", "unit")
        self.tree = ttk.Treeview(tree_frame, columns=columns, show='headings', selectmode="extended", height=10)
        self.tree.heading("category", text="Category", anchor=tk.W); self.tree.column("category", width=180, stretch=tk.NO, anchor=tk.W)
        self.tree.heading("parameter", text="Parameter", anchor=tk.W); self.tree.column("parameter", width=280, stretch=tk.NO, anchor=tk.W)
        self.tree.heading("value", text="Value", anchor=tk.E); self.tree.column("value", width=150, stretch=tk.NO, anchor=tk.E)
        self.tree.heading("unit", text="Unit", anchor=tk.W); self.tree.column("unit", width=100, stretch=tk.NO, anchor=tk.W)
        vsb = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree.yview); hsb = ttk.Scrollbar(tree_frame, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set); vsb.pack(side='right', fill='y'); hsb.pack(side='bottom', fill='x'); self.tree.pack(side='left', fill='both', expand=True)
        self.tree.tag_configure('oddrow', background='#E8E8E8'); self.tree.tag_configure('evenrow', background='#FFFFFF'); self.tree.tag_configure('separator', background='#B0B0B0', foreground='#B0B0B0')

        self.ch_frame = ttk.LabelFrame(self.right_pane, text="Curtin-Hammett Pathway Analysis", padding=10)
        self.ch_frame.pack(fill=tk.X, pady=10, side=tk.BOTTOM)
        
        # Row 0: Reactant, Product, TS Listbox Label
        ttk.Label(self.ch_frame, text="Reactant:").grid(row=0, column=0, padx=5, pady=3, sticky='w')
        self.ch_reactant_combo = ttk.Combobox(self.ch_frame, textvariable=self.ch_reactant_var, state="readonly", width=30)
        self.ch_reactant_combo.grid(row=0, column=1, padx=5, pady=3, sticky='ew')
        self.ch_reactant_combo.bind('<<ComboboxSelected>>', lambda e: self.check_ch_button_state())
        
        ttk.Label(self.ch_frame, text="Transition State(s):").grid(row=0, column=2, padx=(10,2), pady=3, sticky='nw')
        self.ch_ts_listbox_frame = ttk.Frame(self.ch_frame); 
        self.ch_ts_listbox_frame.grid(row=0, column=3, padx=5, pady=3, sticky='ewns', rowspan=2) # Span 2 rows for listbox
        self.ch_ts_listbox = tk.Listbox(self.ch_ts_listbox_frame, selectmode=tk.EXTENDED, exportselection=False, height=3, width=30)
        self.ch_ts_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True); 
        ts_list_sb = ttk.Scrollbar(self.ch_ts_listbox_frame, orient="vertical", command=self.ch_ts_listbox.yview)
        ts_list_sb.pack(side=tk.RIGHT, fill=tk.Y); self.ch_ts_listbox.config(yscrollcommand=ts_list_sb.set)
        self.ch_ts_listbox.bind('<<ListboxSelect>>', lambda e: self.check_ch_button_state())
        
        ch_ts_ctrl_frame = ttk.Frame(self.ch_frame); 
        ch_ts_ctrl_frame.grid(row=0, column=4, padx=2, pady=3, rowspan=2, sticky='ns') # Span 2 rows for buttons
        ttk.Button(ch_ts_ctrl_frame, text="Add TS", command=self.add_ch_ts_files, width=10).pack(fill=tk.X, pady=1)
        ttk.Button(ch_ts_ctrl_frame, text="Remove Sel.", command=self.remove_ch_ts_file, width=10).pack(fill=tk.X, pady=1)
        ttk.Button(ch_ts_ctrl_frame, text="Clear All", command=self.clear_ch_ts_files, width=10).pack(fill=tk.X, pady=1)

        # Row 1: Product
        ttk.Label(self.ch_frame, text="Product (Optional):").grid(row=1, column=0, padx=5, pady=3, sticky='w')
        self.ch_product_combo = ttk.Combobox(self.ch_frame, textvariable=self.ch_product_var, state="readonly", width=30)
        self.ch_product_combo.grid(row=1, column=1, padx=5, pady=3, sticky='ew')
        self.ch_product_combo.bind('<<ComboboxSelected>>', lambda e: self.ch_product_combo.selection_clear() if self.ch_product_var.get() == "(None)" else None)

        # Row 2: CH Plot Options
        ch_plot_options_frame = ttk.Frame(self.ch_frame)
        ch_plot_options_frame.grid(row=2, column=0, columnspan=5, sticky='ew', pady=(5,0))
        
        ttk.Label(ch_plot_options_frame, text="Line Style:").pack(side=tk.LEFT, padx=(0,2))
        self.ch_plot_line_style_combo = ttk.Combobox(ch_plot_options_frame, textvariable=self.ch_plot_line_style_var, values=self.line_styles, state="readonly", width=7)
        self.ch_plot_line_style_combo.pack(side=tk.LEFT, padx=(0,5))
        ttk.Label(ch_plot_options_frame, text="Width:").pack(side=tk.LEFT, padx=(0,2))
        self.ch_plot_line_width_entry = ttk.Entry(ch_plot_options_frame, textvariable=self.ch_plot_line_width_var, width=4)
        self.ch_plot_line_width_entry.pack(side=tk.LEFT, padx=(0,5))
        ttk.Label(ch_plot_options_frame, text="Color:").pack(side=tk.LEFT, padx=(0,2))
        self.ch_plot_line_color_combo = ttk.Combobox(ch_plot_options_frame, textvariable=self.ch_plot_line_color_var, values=self.plot_colors, state="readonly", width=8)
        self.ch_plot_line_color_combo.pack(side=tk.LEFT, padx=(0,10))

        self.enable_ch_axis_break_cb = ttk.Checkbutton(ch_plot_options_frame, text="Break Y-Axis", variable=self.enable_ch_axis_break_var, command=self._toggle_ch_axis_break_entries)
        self.enable_ch_axis_break_cb.pack(side=tk.LEFT, padx=(0,5))
        ttk.Label(ch_plot_options_frame, text="From:").pack(side=tk.LEFT)
        self.ch_axis_break_from_entry = ttk.Entry(ch_plot_options_frame, textvariable=self.ch_axis_break_from_var, width=7, state=tk.DISABLED)
        self.ch_axis_break_from_entry.pack(side=tk.LEFT, padx=(0,5))
        ttk.Label(ch_plot_options_frame, text="To:").pack(side=tk.LEFT)
        self.ch_axis_break_to_entry = ttk.Entry(ch_plot_options_frame, textvariable=self.ch_axis_break_to_var, width=7, state=tk.DISABLED)
        self.ch_axis_break_to_entry.pack(side=tk.LEFT, padx=(0,10))

        # Row 3: Analyze and Plot Buttons
        ch_action_buttons_frame = ttk.Frame(self.ch_frame)
        ch_action_buttons_frame.grid(row=3, column=0, columnspan=5, pady=5, sticky='ew')
        self.ch_analyze_button = ttk.Button(ch_action_buttons_frame, text="Analyze Pathways", command=self.calculate_and_display_ch_params, state=tk.DISABLED)
        self.ch_analyze_button.pack(side=tk.LEFT, padx=(0,5))
        self.ch_plot_button = ttk.Button(ch_action_buttons_frame, text="Plot CH Profile", command=self.plot_curtin_hammett_profile, state=tk.DISABLED)
        self.ch_plot_button.pack(side=tk.LEFT, padx=5)


        # Row 4: Results Tree
        self.ch_results_tree = ttk.Treeview(self.ch_frame, columns=("pathway", "dG_act", "dH_act", "TdS_act", "dS_act", "ratio", "details"), show='headings', height=4) # Adjusted height
        self.ch_results_tree.grid(row=4, column=0, columnspan=5, pady=5, sticky='ewns')
        ch_res_sb_v = ttk.Scrollbar(self.ch_frame, orient="vertical", command=self.ch_results_tree.yview)
        ch_res_sb_v.grid(row=4, column=5, sticky='ns'); self.ch_results_tree.configure(yscrollcommand=ch_res_sb_v.set)
        
        self.ch_results_tree.heading("pathway", text="Pathway / Product via TS", anchor=tk.W); self.ch_results_tree.column("pathway", width=220, anchor=tk.W, stretch=tk.YES)
        self.ch_results_tree.heading("dG_act", text="ΔG‡", anchor=tk.E); self.ch_results_tree.column("dG_act", width=80, anchor=tk.E, stretch=tk.NO)
        self.ch_results_tree.heading("dH_act", text="ΔH‡", anchor=tk.E); self.ch_results_tree.column("dH_act", width=80, anchor=tk.E, stretch=tk.NO)
        self.ch_results_tree.heading("TdS_act", text="-TΔS‡", anchor=tk.E); self.ch_results_tree.column("TdS_act", width=80, anchor=tk.E, stretch=tk.NO)
        self.ch_results_tree.heading("dS_act", text="ΔS‡", anchor=tk.E); self.ch_results_tree.column("dS_act", width=100, anchor=tk.E, stretch=tk.NO)
        self.ch_results_tree.heading("ratio", text="Ratio (%)", anchor=tk.E); self.ch_results_tree.column("ratio", width=80, anchor=tk.E, stretch=tk.NO)
        self.ch_results_tree.heading("details", text="Units/Info", anchor=tk.W); self.ch_results_tree.column("details", width=150, anchor=tk.W, stretch=tk.YES)

        self.ch_frame.grid_columnconfigure(1, weight=1)
        self.ch_frame.grid_columnconfigure(3, weight=1) # Adjust weight for TS listbox column
        self.ch_frame.grid_rowconfigure(4, weight=1) # Allow results tree to expand

    def _toggle_axis_break_entries_tab1(self):
        if self.enable_axis_break_var.get():
            self.axis_break_from_entry.config(state=tk.NORMAL)
            self.axis_break_to_entry.config(state=tk.NORMAL)
        else:
            self.axis_break_from_entry.config(state=tk.DISABLED)
            self.axis_break_to_entry.config(state=tk.DISABLED)
            self.axis_break_from_var.set("")
            self.axis_break_to_var.set("")
            
    def _toggle_axis_break_entries_tab2(self):
        if self.enable_allxyz_axis_break_var.get():
            self.allxyz_axis_break_from_entry.config(state=tk.NORMAL)
            self.allxyz_axis_break_to_entry.config(state=tk.NORMAL)
        else:
            self.allxyz_axis_break_from_entry.config(state=tk.DISABLED)
            self.allxyz_axis_break_to_entry.config(state=tk.DISABLED)
            self.allxyz_axis_break_from_var.set("")
            self.allxyz_axis_break_to_var.set("")
            
    def _toggle_ch_axis_break_entries(self): # New method for CH plot axis break
        if self.enable_ch_axis_break_var.get():
            self.ch_axis_break_from_entry.config(state=tk.NORMAL)
            self.ch_axis_break_to_entry.config(state=tk.NORMAL)
        else:
            self.ch_axis_break_from_entry.config(state=tk.DISABLED)
            self.ch_axis_break_to_entry.config(state=tk.DISABLED)
            self.ch_axis_break_from_var.set("")
            self.ch_axis_break_to_var.set("")


    def setup_tab2_ui(self):
        master = self.tab2
        load_allxyz_frame = ttk.Frame(master, padding=10); load_allxyz_frame.pack(fill=tk.X)
        tk.Button(load_allxyz_frame, text="Load .allxyz / multi-frame .xyz File(s)", command=self.load_allxyz_files).pack(side=tk.LEFT, padx=5)
        self.allxyz_filepath_label = tk.Label(load_allxyz_frame, text="No .allxyz file loaded", relief=tk.SUNKEN, width=60, anchor='w')
        self.allxyz_filepath_label.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)

        allxyz_paned_window = tk.PanedWindow(master, orient=tk.HORIZONTAL, sashrelief=tk.RAISED) 
        allxyz_paned_window.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        allxyz_list_frame = ttk.Frame(allxyz_paned_window); allxyz_paned_window.add(allxyz_list_frame, width=250, minsize=200)
        ttk.Label(allxyz_list_frame, text="Loaded .allxyz Files:").pack(anchor='w')
        self.allxyz_listbox = tk.Listbox(allxyz_list_frame, selectmode=tk.SINGLE, exportselection=False, height=10)
        self.allxyz_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        allxyz_list_sb = ttk.Scrollbar(allxyz_list_frame, orient="vertical", command=self.allxyz_listbox.yview)
        allxyz_list_sb.pack(side=tk.RIGHT, fill=tk.Y); self.allxyz_listbox.config(yscrollcommand=allxyz_list_sb.set)
        self.allxyz_listbox.bind("<<ListboxSelect>>", self.on_allxyz_file_select)

        allxyz_details_frame = ttk.Frame(allxyz_paned_window); allxyz_paned_window.add(allxyz_details_frame)
        
        allxyz_plot_options_frame = ttk.LabelFrame(allxyz_details_frame, text="Plot Options (.allxyz Profile)", padding=5)
        allxyz_plot_options_frame.pack(fill=tk.X, pady=5)

        # Plot Style options for allxyz plot
        style_controls_frame = ttk.Frame(allxyz_plot_options_frame)
        style_controls_frame.pack(fill=tk.X, pady=2)
        ttk.Label(style_controls_frame, text="Line Style:").grid(row=0, column=0, padx=(0,2), pady=2, sticky=tk.W)
        self.allxyz_plot_line_style_combo = ttk.Combobox(style_controls_frame, textvariable=self.allxyz_plot_line_style_var, values=self.line_styles, state="readonly", width=7)
        self.allxyz_plot_line_style_combo.grid(row=0, column=1, padx=(0,5), pady=2, sticky=tk.W)
        ttk.Label(style_controls_frame, text="Width:").grid(row=0, column=2, padx=(5,2), pady=2, sticky=tk.W)
        self.allxyz_plot_line_width_entry = ttk.Entry(style_controls_frame, textvariable=self.allxyz_plot_line_width_var, width=4)
        self.allxyz_plot_line_width_entry.grid(row=0, column=3, padx=(0,5), pady=2, sticky=tk.W)
        ttk.Label(style_controls_frame, text="Color:").grid(row=0, column=4, padx=(5,2), pady=2, sticky=tk.W)
        self.allxyz_plot_line_color_combo = ttk.Combobox(style_controls_frame, textvariable=self.allxyz_plot_line_color_var, values=self.plot_colors, state="readonly", width=8)
        self.allxyz_plot_line_color_combo.grid(row=0, column=5, padx=(0,5), pady=2, sticky=tk.W)
        
        # Axis Break options for allxyz plot
        allxyz_axis_break_frame = ttk.Frame(allxyz_plot_options_frame)
        allxyz_axis_break_frame.pack(fill=tk.X, pady=(5,0))
        self.enable_allxyz_axis_break_cb = ttk.Checkbutton(allxyz_axis_break_frame, text="Enable Y-Axis Break", variable=self.enable_allxyz_axis_break_var, command=self._toggle_axis_break_entries_tab2)
        self.enable_allxyz_axis_break_cb.grid(row=0, column=0, padx=(0,5), pady=2, sticky=tk.W)
        ttk.Label(allxyz_axis_break_frame, text="Break From (kcal/mol):").grid(row=0, column=1, padx=(0,2), pady=2, sticky=tk.W)
        self.allxyz_axis_break_from_entry = ttk.Entry(allxyz_axis_break_frame, textvariable=self.allxyz_axis_break_from_var, width=8, state=tk.DISABLED)
        self.allxyz_axis_break_from_entry.grid(row=0, column=2, padx=(0,5), pady=2, sticky=tk.W)
        ttk.Label(allxyz_axis_break_frame, text="Break To (kcal/mol):").grid(row=0, column=3, padx=(0,2), pady=2, sticky=tk.W)
        self.allxyz_axis_break_to_entry = ttk.Entry(allxyz_axis_break_frame, textvariable=self.allxyz_axis_break_to_var, width=8, state=tk.DISABLED)
        self.allxyz_axis_break_to_entry.grid(row=0, column=4, padx=(0,0), pady=2, sticky=tk.W)
        
        self.highest_e_info_label = ttk.Label(allxyz_details_frame, text="Highest Energy Frame: N/A", wraplength=500)
        self.highest_e_info_label.pack(pady=5, anchor='w', before=allxyz_plot_options_frame) # Ensure it's above plot options

        allxyz_action_frame = ttk.Frame(allxyz_details_frame); allxyz_action_frame.pack(fill=tk.X, pady=5)
        self.save_highest_e_button = ttk.Button(allxyz_action_frame, text="Save Highest E Geometry (.xyz)", command=self.save_highest_energy_geometry, state=tk.DISABLED)
        self.save_highest_e_button.pack(side=tk.LEFT, padx=5)
        self.plot_allxyz_button = ttk.Button(allxyz_action_frame, text="Plot Profile for this .allxyz", command=self.plot_allxyz_profile, state=tk.DISABLED)
        self.plot_allxyz_button.pack(side=tk.LEFT, padx=5)
        self.copy_allxyz_table_button = ttk.Button(allxyz_action_frame, text="Copy Frame Table (TSV)", command=self.copy_allxyz_table_to_clipboard)
        self.copy_allxyz_table_button.pack(side=tk.LEFT, padx=5)


        allxyz_tree_frame = ttk.Frame(allxyz_details_frame); allxyz_tree_frame.pack(fill=tk.BOTH, expand=True, pady=5)
        self.allxyz_details_tree = ttk.Treeview(allxyz_tree_frame, columns=("frame", "energy_eh", "rel_energy_kcal"), show='headings', height=15)
        self.allxyz_details_tree.heading("frame", text="Frame #", anchor=tk.W); self.allxyz_details_tree.column("frame", width=80, anchor=tk.W, stretch=tk.NO)
        self.allxyz_details_tree.heading("energy_eh", text="Energy (Eh)", anchor=tk.E); self.allxyz_details_tree.column("energy_eh", width=180, anchor=tk.E, stretch=tk.NO)
        self.allxyz_details_tree.heading("rel_energy_kcal", text="Rel. E (kcal/mol)", anchor=tk.E); self.allxyz_details_tree.column("rel_energy_kcal", width=180, anchor=tk.E, stretch=tk.NO)
        allxyz_tree_vsb = ttk.Scrollbar(allxyz_tree_frame, orient="vertical", command=self.allxyz_details_tree.yview)
        allxyz_tree_hsb = ttk.Scrollbar(allxyz_tree_frame, orient="horizontal", command=self.allxyz_details_tree.xview)
        self.allxyz_details_tree.configure(yscrollcommand=allxyz_tree_vsb.set, xscrollcommand=allxyz_tree_hsb.set)
        allxyz_tree_vsb.pack(side='right', fill='y'); allxyz_tree_hsb.pack(side='bottom', fill='x'); self.allxyz_details_tree.pack(fill=tk.BOTH, expand=True)
        self.allxyz_details_tree.tag_configure('highest_e', background='lightyellow')

    # --- Methods for TAB 1 (Thermo Analysis & Curtin-Hammett) ---
    def apply_treeview_row_tags(self, treeview_widget): # Made generic
        children = treeview_widget.get_children(''); data_row_index = 0
        for item_id in children:
            cv = treeview_widget.item(item_id, 'values')
            if cv and (cv[0] == "---" or cv[0] in ["Error", "Info"]):
                if cv[0] == "---": treeview_widget.item(item_id, tags=('separator',))
            else:
                treeview_widget.item(item_id, tags=(self.row_tags[data_row_index % 2],))
                data_row_index += 1

    def update_plot_reference_selector(self):
        display_names = [os.path.basename(path) for path in self.loaded_file_paths if path in self.loaded_files_data and self.loaded_files_data[path] and \
                        (self.loaded_files_data[path].get("final_gibbs_free_energy_g_eh") or self.loaded_files_data[path].get("electronic_energy_eh"))]
        plot_ref_options = ["[Use Overall Lowest Energy]"] + sorted(list(set(display_names)))
        self.plot_reference_combo['values'] = plot_ref_options
        current_ref = self.plot_reference_var.get()
        if not current_ref or current_ref not in plot_ref_options : self.plot_reference_var.set(plot_ref_options[0])

    def check_plot_button_state(self): 
        # Check if there's any plottable energy for the selected type across all loaded files
        plottable_found = False
        selected_plot_type_key_base = self.plot_energy_type_var.get()
        for data in self.loaded_files_data.values():
            if data:
                energy_key = self._get_energy_key_for_plot(selected_plot_type_key_base, data)
                if energy_key == "_e_plus_zpe_calculated":
                    if data.get("electronic_energy_eh") is not None and data.get("zero_point_energy_eh") is not None:
                        plottable_found = True; break
                    elif data.get("electronic_energy_eh") is not None: # Fallback for E+ZPE if ZPE missing
                        plottable_found = True; break
                elif data.get(energy_key) is not None:
                    plottable_found = True; break
                elif data.get("electronic_energy_eh") is not None: # Ultimate fallback to E_el
                    plottable_found = True; break
        
        if plottable_found and MATPLOTLIB_AVAILABLE:
            self.plot_energy_button.config(state=tk.NORMAL)
        else:
            self.plot_energy_button.config(state=tk.DISABLED)

    def load_files(self):
        filepaths_tuple = filedialog.askopenfilenames(
            title="Select ORCA Output File(s)",
            filetypes=(("ORCA files", "*.out *.log *.xyz"), ("All files", "*.*"))
        )
        if not filepaths_tuple:
            return

        # Store current selection to try and restore it
        old_selected_path = self.current_selected_filepath 
        
        newly_added_to_main_list_paths = []

        for filepath in filepaths_tuple:
            if filepath not in self.loaded_file_paths:
                self.loaded_file_paths.append(filepath)
                newly_added_to_main_list_paths.append(filepath)
            
            # Always parse (or re-parse) and store/update the data
            parsed_data = parse_orca_thermo_block(filepath)
            if parsed_data:
                self.loaded_files_data[filepath] = parsed_data
            else:
                # If parsing failed, ensure no stale data exists for this path
                if filepath in self.loaded_files_data:
                    del self.loaded_files_data[filepath]
                print(f"Warning: Parsing failed for {os.path.basename(filepath)}. It will be marked in the list.")

        # Re-populate the listbox from the (potentially modified/extended) self.loaded_file_paths
        self.file_listbox.delete(0, tk.END)
        for f_path in self.loaded_file_paths:
            filename = os.path.basename(f_path)
            self.file_listbox.insert(tk.END, filename)
            # Color red if parsing failed or data is None for this path
            if f_path not in self.loaded_files_data or self.loaded_files_data[f_path] is None:
                try:
                    # Find index of filename in listbox items to color it
                    # This assumes filenames are unique in the listbox, which they should be if paths are unique
                    idx_to_color = self.file_listbox.get(0, tk.END).index(filename)
                    self.file_listbox.itemconfig(idx_to_color, {'fg': 'red'})
                except ValueError:
                    # This case should ideally not be reached if f_path is from self.loaded_file_paths
                    # and was just inserted into the listbox.
                    print(f"Debug: Could not find {filename} in listbox to color red.")
                    pass 

        # Update UI elements that depend on the full list of files
        self.update_plot_reference_selector()
        self.update_ch_selectors() 
        self.check_plot_button_state()
        self.check_ch_button_state() 

        # Reset the detailed display for a single file; selection logic below will trigger its update
        self.filepath_label.config(text="Select a file for display")
        for item in self.tree.get_children(): self.tree.delete(item)
        self.display_button.config(state=tk.DISABLED)

        # Logic for setting selection in the listbox
        if self.file_listbox.size() > 0:
            idx_to_select = -1
            if newly_added_to_main_list_paths: # Prioritize selecting the first *newly* added file
                try:
                    # Find index in self.loaded_file_paths, which matches listbox order
                    idx_to_select = self.loaded_file_paths.index(newly_added_to_main_list_paths[0])
                except ValueError:
                    pass 
            elif old_selected_path and old_selected_path in self.loaded_file_paths: # Try to re-select old one
                try:
                    idx_to_select = self.loaded_file_paths.index(old_selected_path)
                except ValueError:
                    pass
            
            if idx_to_select != -1 and idx_to_select < self.file_listbox.size():
                self.file_listbox.selection_set(idx_to_select)
                self.file_listbox.see(idx_to_select) 
                self.file_listbox.activate(idx_to_select)
            else: # Default to first item if no specific selection could be made
                self.file_listbox.selection_set(0)
                self.file_listbox.see(0)
                self.file_listbox.activate(0)
            
            self.file_listbox.event_generate("<<ListboxSelect>>") 
        else:
            # No files in listbox, clear displays properly
            self.current_selected_filepath = None
            self.filepath_label.config(text="No file selected for display")
            self.display_button.config(state=tk.DISABLED)
            for item in self.tree.get_children(): self.tree.delete(item)
            self.tree.insert('', tk.END, values=("Info", "Load ORCA file(s).", "", ""))


    def on_file_select_from_listbox(self, event=None):
        widget = event.widget if event else self.file_listbox; selected_indices = widget.curselection()
        if not selected_indices: 
            # If nothing is selected (e.g. list is empty after removing all files)
            self.current_selected_filepath = None
            self.filepath_label.config(text="No file selected for display")
            for item in self.tree.get_children(): self.tree.delete(item)
            self.tree.insert('', tk.END, values=("Info", "Load ORCA file(s) or select one from the list.", "", ""))
            self.display_button.config(state=tk.DISABLED)
            return

        try: 
            # Ensure index is valid for self.loaded_file_paths
            selected_listbox_index = selected_indices[0]
            if 0 <= selected_listbox_index < len(self.loaded_file_paths):
                 self.current_selected_filepath = self.loaded_file_paths[selected_listbox_index]
            else: # Should not happen if listbox and loaded_file_paths are in sync
                self.filepath_label.config(text="Selection Error."); self.display_button.config(state=tk.DISABLED); return

        except IndexError: 
            self.filepath_label.config(text="Error selecting file."); self.display_button.config(state=tk.DISABLED); return
        
        self.filepath_label.config(text=f"Displaying: {os.path.basename(self.current_selected_filepath)}")
        
        # Check if data exists and is not None (parsing was successful)
        if self.current_selected_filepath in self.loaded_files_data and self.loaded_files_data[self.current_selected_filepath] is not None: 
            self.display_button.config(state=tk.NORMAL); self.refresh_display()
        else:
            for item in self.tree.get_children(): self.tree.delete(item)
            self.tree.insert('', tk.END, values=("Error", f"Data not available or failed to parse for {os.path.basename(self.current_selected_filepath)}.", "", ""))
            self.display_button.config(state=tk.DISABLED)

    def on_selection_change(self): # Covers units and sections for Tab 1 display
        if self.current_selected_filepath and self.current_selected_filepath in self.loaded_files_data:
            self.refresh_display()
        self.check_ch_button_state()

    def trigger_refresh_display(self):
        if not self.current_selected_filepath: messagebox.showwarning("No File Selected", "Select a file from list."); return
        self.refresh_display()

    def refresh_display(self):
        if self.current_selected_filepath and self.current_selected_filepath in self.loaded_files_data and self.loaded_files_data[self.current_selected_filepath] is not None:
            populate_treeview_widget(self.tree, self.loaded_files_data[self.current_selected_filepath], self)
        else:
            for item in self.tree.get_children(): self.tree.delete(item)
            if self.current_selected_filepath: self.tree.insert('', tk.END, values=("Info", f"No data for {os.path.basename(self.current_selected_filepath)} (parsing may have failed).", "", ""))
            elif self.file_listbox.size() > 0: self.tree.insert('', tk.END, values=("Info", "Select a file to view.", "", ""))
            else: self.tree.insert('', tk.END, values=("Info", "Load ORCA file(s).", "", ""))
        self.check_ch_button_state()
        self.check_plot_button_state() 

    def copy_thermo_table_to_clipboard(self):
        if not self.tree.get_children() or \
           (len(self.tree.get_children()) == 1 and self.tree.item(self.tree.get_children()[0], 'values')[0] in ["Error", "Info"]):
            messagebox.showinfo("Clipboard", "Thermo table empty."); return
        header = ["Category", "Parameter", "Value", "Unit"]; tsv_data = ["\t".join(header)]
        for item_id in self.tree.get_children():
            values = self.tree.item(item_id, 'values')
            if values and values[0] == "---": continue
            tsv_data.append("\t".join([str(v) if v is not None else "" for v in values]))
        full_tsv_string = "\n".join(tsv_data)
        try: self.master.clipboard_clear(); self.master.clipboard_append(full_tsv_string); messagebox.showinfo("Clipboard", "Thermo table copied (TSV).")
        except tk.TclError: messagebox.showerror("Clipboard Error", "Could not access clipboard.")

    def _plot_profile_matplotlib(self, x_coords, y_values_kcal, x_labels, y_label_main, title_text, 
                                 ref_label_for_plot, types_plotted_for_annotations, 
                                 enable_break, break_from_str, break_to_str, 
                                 line_style, line_width_str, line_color, # New style parameters
                                 is_allxyz_plot=False):
        """Helper function to handle Matplotlib plotting with optional axis break and styles."""
        if not MATPLOTLIB_AVAILABLE:
            messagebox.showerror("Plot Error", "Matplotlib library is not installed. Please install it to use plotting features.")
            return

        fig = None
        ax_main = None 
        ax_top = None  

        try:
            line_width = float(line_width_str)
        except ValueError:
            messagebox.showerror("Plot Error", f"Invalid line width: '{line_width_str}'. Must be a number.", parent=self.master)
            return

        try:
            break_from_val = float(break_from_str) if enable_break and break_from_str else None
            break_to_val = float(break_to_str) if enable_break and break_to_str else None
        except ValueError:
            messagebox.showerror("Plot Error", "Invalid axis break values. Please enter numbers.", parent=self.master) 
            return

        if enable_break and break_from_val is not None and break_to_val is not None and break_from_val < break_to_val:
            fig, (ax_top, ax_main) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [1, 1]}) 
            fig.subplots_adjust(hspace=0.1)  

            ax_main.set_ylim(min(y_values_kcal) - 5, break_from_val)  
            ax_top.set_ylim(break_to_val, max(y_values_kcal) + 5)    

            ax_main.spines['top'].set_visible(False)
            ax_top.spines['bottom'].set_visible(False)
            ax_top.tick_params(axis='x', which='both', bottom=False) 

            d = .015  
            kwargs_break = dict(transform=ax_main.transAxes, color='k', clip_on=False, linewidth=1)
            ax_main.plot((-d, +d), (1 - d, 1 + d), **kwargs_break)        
            ax_main.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs_break)  

            kwargs_break.update(transform=ax_top.transAxes)
            ax_top.plot((-d, +d), (0 - d, 0 + d), **kwargs_break)      
            ax_top.plot((1 - d, 1 + d), (0 - d, 0 + d), **kwargs_break) 
            
            marker_s = 25 if is_allxyz_plot else 20
            marker_ew = 2 if is_allxyz_plot else 3
            
            ax_main.plot(x_coords, y_values_kcal, marker='_', linestyle=line_style, color=line_color, markersize=marker_s, markeredgewidth=marker_ew, linewidth=line_width)
            ax_top.plot(x_coords, y_values_kcal, marker='_', linestyle=line_style, color=line_color, markersize=marker_s, markeredgewidth=marker_ew, linewidth=line_width)

            fig.text(0.02, 0.5, f"Relative {y_label_main} (kcal/mol)\n(Zeroed at: {ref_label_for_plot})", va='center', rotation='vertical', fontsize=10)
            ax_main.set_xlabel("Reaction Coordinate" if not is_allxyz_plot else "Frame Number", fontsize=10)
            
            for i, txt_val in enumerate(y_values_kcal):
                label_text = f"{txt_val:.2f}"
                if types_plotted_for_annotations and i < len(types_plotted_for_annotations):
                    label_text += f"\n({types_plotted_for_annotations[i]})"

                if y_values_kcal[i] <= break_from_val: 
                    ax_main.annotate(label_text, (x_coords[i], y_values_kcal[i]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=8)
                elif y_values_kcal[i] >= break_to_val: 
                     ax_top.annotate(label_text, (x_coords[i], y_values_kcal[i]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=8)

            ax_main.grid(axis='y', linestyle=':', alpha=0.8, color='gray')
            ax_top.grid(axis='y', linestyle=':', alpha=0.8, color='gray')
            ax_main.axhline(0, color='black', linewidth=0.5, linestyle='--') 
            if 0 > break_to_val : 
                 ax_top.axhline(0, color='black', linewidth=0.5, linestyle='--')


        else: # No axis break
            fig, ax_main = plt.subplots(figsize=(max(7, len(x_labels) * 0.9), 6.5))
            marker_s = 25 if is_allxyz_plot else 20
            marker_ew = 2 if is_allxyz_plot else 3

            ax_main.plot(x_coords, y_values_kcal, marker='_', linestyle=line_style, color=line_color, markersize=marker_s, markeredgewidth=marker_ew, linewidth=line_width)
            ax_main.set_ylabel(f"Relative {y_label_main} (kcal/mol)\n(Zeroed at: {ref_label_for_plot})", fontsize=10)
            ax_main.set_xlabel("Reaction Coordinate" if not is_allxyz_plot else "Frame Number", fontsize=10)
            ax_main.grid(axis='y', linestyle=':', alpha=0.8, color='gray')
            ax_main.axhline(0, color='black', linewidth=0.5, linestyle='--')

            for i, txt_val in enumerate(y_values_kcal):
                label_text = f"{txt_val:.2f}"
                if types_plotted_for_annotations and i < len(types_plotted_for_annotations):
                     label_text += f"\n({types_plotted_for_annotations[i]})"
                ax_main.annotate(label_text, (x_coords[i], y_values_kcal[i]), textcoords="offset points", xytext=(0,10), ha='center', fontsize=8)

        if ax_main is not None: 
            ax_main.set_xticks(x_coords)
            ax_main.set_xticklabels(x_labels, rotation=40, ha="right", fontsize=9)
            if len(x_labels) > 20: 
                ax_main.xaxis.set_major_locator(MaxNLocator(integer=True, prune='both', nbins=15))


        plt.suptitle(title_text, fontsize=12, fontweight='bold')
        plt.tight_layout(rect=[0.03, 0.03, 1, 0.95]) 
        plt.show()


    def plot_energy_profile(self):
        if not MATPLOTLIB_AVAILABLE: messagebox.showerror("Plot Error", "Matplotlib not installed."); return
        plot_data = []
        requested_type_str = self.plot_energy_type_var.get()
        y_axis_label_main = requested_type_str.split('(')[0].strip()

        # Get plot style parameters
        line_style = self.plot_line_style_var.get()
        line_width_str = self.plot_line_width_var.get()
        line_color = self.plot_line_color_var.get()


        for path in self.loaded_file_paths: 
            if path in self.loaded_files_data and self.loaded_files_data[path]:
                data = self.loaded_files_data[path]; energy_to_plot = None; type_plotted_label = "N/A"
                
                if requested_type_str == "Gibbs Free Energy (G)": 
                    energy_to_plot = data.get("final_gibbs_free_energy_g_eh"); type_plotted_label = "G"
                elif requested_type_str == "Enthalpy (H)": 
                    energy_to_plot = data.get("total_enthalpy_h_eh"); type_plotted_label = "H"
                elif requested_type_str == "E_el + ZPVE":
                    e_el = data.get("electronic_energy_eh"); zpve = data.get("zero_point_energy_eh")
                    if e_el is not None and zpve is not None: energy_to_plot = e_el + zpve; type_plotted_label = "E+ZPE"
                    elif e_el is not None: energy_to_plot = e_el; type_plotted_label = "E_el (ZPE miss)"
                elif requested_type_str == "Electronic Energy (E_el)": 
                    energy_to_plot = data.get("electronic_energy_eh"); type_plotted_label = "E_el"
                
                if energy_to_plot is None and data.get("electronic_energy_eh") is not None: 
                    energy_to_plot = data.get("electronic_energy_eh"); type_plotted_label = f"E_el (fallbk for {type_plotted_label if type_plotted_label != 'N/A' else requested_type_str.split('(')[0].strip()})"
                
                if energy_to_plot is not None:
                    plot_data.append({"name": os.path.basename(path).rsplit('.',1)[0], "energy_h": energy_to_plot, "full_path": path, "type_plotted": type_plotted_label})
        
        if len(plot_data) < 1: messagebox.showinfo("Energy Profile", "Not enough data for selected type."); return

        names = [item["name"] for item in plot_data]; energies_hartree = [item["energy_h"] for item in plot_data]
        types_plotted = [item["type_plotted"] for item in plot_data]
        selected_ref_display_name = self.plot_reference_var.get(); ref_energy_hartree = None; ref_label_for_plot = "Overall Lowest"

        if selected_ref_display_name == "[Use Overall Lowest Energy]" or not selected_ref_display_name:
            if not energies_hartree: messagebox.showerror("Plot Error", "No energies."); return
            ref_energy_hartree = min(energies_hartree); min_idx = energies_hartree.index(ref_energy_hartree)
            ref_label_for_plot = f"Lowest: {names[min_idx]} ({types_plotted[min_idx]}={ref_energy_hartree:.6f} Ha)"
        else:
            ref_item_data = next((item for item in plot_data if item["name"] == selected_ref_display_name.rsplit('.',1)[0]), None) 
            if ref_item_data: ref_energy_hartree = ref_item_data["energy_h"]; ref_type_plotted_for_ref = ref_item_data["type_plotted"]; ref_label_for_plot = f"{selected_ref_display_name} ({ref_type_plotted_for_ref}={ref_energy_hartree:.6f} Ha)"
            else:
                if not energies_hartree: messagebox.showerror("Plot Error", "No energies."); return
                messagebox.showwarning("Plot Ref Warning", f"Ref '{selected_ref_display_name}' data problem. Defaulting to lowest.")
                ref_energy_hartree = min(energies_hartree); min_idx = energies_hartree.index(ref_energy_hartree); ref_label_for_plot = f"Lowest: {names[min_idx]} ({types_plotted[min_idx]}={ref_energy_hartree:.6f} Ha)"
        if ref_energy_hartree is None: messagebox.showerror("Plot Error", "Cannot determine reference energy."); return
        relative_energies_kcal = [(e - ref_energy_hartree) * HARTREE_TO_KCAL_PER_MOL for e in energies_hartree]

        self._plot_profile_matplotlib(
            x_coords=range(len(names)),
            y_values_kcal=relative_energies_kcal,
            x_labels=names,
            y_label_main=y_axis_label_main,
            title_text=f"Energy Profile ({y_axis_label_main})",
            ref_label_for_plot=ref_label_for_plot,
            types_plotted_for_annotations=types_plotted,
            enable_break=self.enable_axis_break_var.get(),
            break_from_str=self.axis_break_from_var.get(),
            break_to_str=self.axis_break_to_var.get(),
            line_style=line_style,
            line_width_str=line_width_str,
            line_color=line_color,
            is_allxyz_plot=False
        )


    # --- Curtin-Hammett Specific Methods ---
    def _get_ch_filepath_from_display_name(self, display_name):
        if not display_name or display_name == "(None)": return None
        return next((path for path in self.loaded_file_paths if os.path.basename(path) == display_name), None)

    def add_ch_ts_files(self):
        filepaths_tuple = filedialog.askopenfilenames(title="Add Transition State File(s) for Curtin-Hammett",filetypes=(("ORCA files","*.out *.log *.xyz"),("All files","*.*")))
        if not filepaths_tuple: return
        newly_added_to_listbox = False
        for filepath in filepaths_tuple:
            filename = os.path.basename(filepath)
            if filepath not in self.loaded_files_data or self.loaded_files_data[filepath] is None: # Check if data is missing or parsing failed
                parsed_data = parse_orca_thermo_block(filepath)
                if parsed_data:
                    if filepath not in self.loaded_file_paths: 
                        self.loaded_file_paths.append(filepath) # Add to main list if new
                        # Update main listbox if it was a completely new file not previously loaded
                        if filename not in self.file_listbox.get(0, tk.END):
                           self.file_listbox.insert(tk.END, filename) # Add to main GUI list
                           # If it's a new file, color it red if parsing failed (though parsed_data should be True here)
                           if not parsed_data: # Should not happen if parsed_data is True
                               try:
                                   idx = self.file_listbox.get(0, tk.END).index(filename)
                                   self.file_listbox.itemconfig(idx, {'fg': 'red'})
                               except ValueError: pass


                    self.loaded_files_data[filepath] = parsed_data; 
                    self.update_ch_selectors() # Ensure CH selectors are aware of new data
                else: messagebox.showwarning("TS Parse Warning", f"Could not parse TS file: {filename}. It will not be available for CH analysis."); continue
            
            # Add to CH TS listbox if not already there AND if it's a valid file (has data)
            if filepath in self.loaded_files_data and self.loaded_files_data[filepath] is not None:
                if filename not in self.ch_ts_listbox.get(0, tk.END): 
                    self.ch_ts_listbox.insert(tk.END, filename); newly_added_to_listbox = True
            else:
                 messagebox.showwarning("TS Add Warning", f"File {filename} was not parsed correctly or is missing data, cannot add to CH TS list.")

        if newly_added_to_listbox: self.check_ch_button_state()
        self.update_plot_reference_selector() # Update plot ref in case new files were added to main list

    def remove_ch_ts_file(self):
        selected_indices = self.ch_ts_listbox.curselection()
        if not selected_indices: return
        for i in sorted(selected_indices, reverse=True): self.ch_ts_listbox.delete(i)
        self.check_ch_button_state()

    def clear_ch_ts_files(self):
        self.ch_ts_listbox.delete(0, tk.END); self.check_ch_button_state()

    def update_ch_selectors(self):
        # Filter for files that have *some* plottable energy, indicating successful basic parsing
        valid_display_names = []
        for path, data in self.loaded_files_data.items():
            if data: # Check if data is not None (i.e., parsing didn't return None)
                 # Check for any of the key energy terms
                if data.get("final_gibbs_free_energy_g_eh") is not None or \
                   data.get("total_enthalpy_h_eh") is not None or \
                   data.get("electronic_energy_eh") is not None:
                    valid_display_names.append(os.path.basename(path))
        valid_display_names.sort()
        
        current_r = self.ch_reactant_var.get(); current_p = self.ch_product_var.get()
        self.ch_reactant_combo['values'] = valid_display_names if valid_display_names else [""]
        self.ch_product_combo['values'] = ["(None)"] + valid_display_names
        
        if current_r in valid_display_names: self.ch_reactant_var.set(current_r)
        elif valid_display_names: self.ch_reactant_var.set(valid_display_names[0])
        else: self.ch_reactant_var.set("")
        
        if current_p == "(None)" or current_p in valid_display_names : self.ch_product_var.set(current_p)
        else: self.ch_product_var.set("(None)")
        
        current_ts_listbox_items = list(self.ch_ts_listbox.get(0, tk.END))
        self.ch_ts_listbox.delete(0, tk.END) # Clear and re-populate to ensure only valid ones remain
        for item_name in current_ts_listbox_items:
            # Ensure the TS file is still in the main loaded list and has valid data
            ts_path_check = self._get_ch_filepath_from_display_name(item_name)
            if ts_path_check and ts_path_check in self.loaded_files_data and self.loaded_files_data[ts_path_check]:
                self.ch_ts_listbox.insert(tk.END, item_name)
        self.check_ch_button_state()

    def check_ch_button_state(self):
        r_selected = self.ch_reactant_var.get() and self.ch_reactant_var.get() != "(None)"
        ts_selected_count = len(self.ch_ts_listbox.get(0, tk.END))
        self.ch_analyze_button.config(state=tk.NORMAL if r_selected and ts_selected_count > 0 else tk.DISABLED)
        # Also update the CH plot button state
        if hasattr(self, 'ch_plot_button'): # Check if button exists (UI built)
            self.ch_plot_button.config(state=tk.NORMAL if r_selected and ts_selected_count > 0 and MATPLOTLIB_AVAILABLE else tk.DISABLED)


    def calculate_and_display_ch_params(self):
        for item in self.ch_results_tree.get_children(): self.ch_results_tree.delete(item)
        r_name_display = self.ch_reactant_var.get(); p_name_display = self.ch_product_var.get()
        ch_ts_display_names = list(self.ch_ts_listbox.get(0, tk.END))
        r_path = self._get_ch_filepath_from_display_name(r_name_display)
        p_path = self._get_ch_filepath_from_display_name(p_name_display) if p_name_display != "(None)" else None
        ts_details_for_calc = [{'path': self._get_ch_filepath_from_display_name(name_disp), 'name': name_disp, 'data': self.loaded_files_data.get(self._get_ch_filepath_from_display_name(name_disp))} for name_disp in ch_ts_display_names]
        ts_details_for_calc = [ts for ts in ts_details_for_calc if ts['path'] and ts['data']]

        if not r_path or r_path not in self.loaded_files_data or self.loaded_files_data[r_path] is None: 
            messagebox.showerror("CH Error", "Reactant data missing or failed to parse."); return
        r_data = self.loaded_files_data[r_path]; r_e, r_e_type, r_temp = self._get_relevant_energy(r_data)
        r_s_au, _, _ = self._get_s_and_ts(r_data) # r_temp from _get_relevant_energy is preferred
        if r_e is None: messagebox.showerror("CH Error", f"Energy for reactant '{r_name_display}' not found."); return
        if r_temp is None: r_temp = 298.15; messagebox.showwarning("CH Warning", "Reactant temp not found, using 298.15K for ratios.")


        def insert_ch_results(pathway, dG_kcal, dH_kcal, minus_TdS_kcal, dS_cal, ratio_percent_str, details_str):
            self.ch_results_tree.insert('',tk.END, values=(pathway,
                f"{dG_kcal:.2f}" if isinstance(dG_kcal, float) else dG_kcal,
                f"{dH_kcal:.2f}" if isinstance(dH_kcal, float) else dH_kcal,
                f"{minus_TdS_kcal:.2f}" if isinstance(minus_TdS_kcal, float) else minus_TdS_kcal,
                f"{dS_cal:.2f}" if isinstance(dS_cal, float) else dS_cal,
                ratio_percent_str, details_str))

        if p_path and p_path in self.loaded_files_data and self.loaded_files_data[p_path] is not None:
            p_data = self.loaded_files_data[p_path]
            p_G, pG_type, p_temp_eff = self._get_relevant_energy(p_data, ['G'])
            p_H, pH_type, _ = self._get_relevant_energy(p_data, ['H','E_el+ZPE','E_el'])
            p_S_au, _, _ = self._get_s_and_ts(p_data)
            r_G_for_rxn, rG_type, _ = self._get_relevant_energy(r_data, ['G'])
            r_H_for_rxn, rH_type, _ = self._get_relevant_energy(r_data, ['H','E_el+ZPE','E_el'])
            dG_rxn, dH_rxn, dS_rxn, TdS_rxn, notes_rxn = "N/A", "N/A", "N/A", "N/A", "Missing data"
            temp_for_rxn = p_temp_eff if p_temp_eff is not None else r_temp # Prioritize product's temp if available

            if p_G is not None and r_G_for_rxn is not None:
                dG_rxn = (p_G - r_G_for_rxn) * HARTREE_TO_KCAL_PER_MOL; notes_rxn = f"ΔG(G_P-G_R)"
                if p_S_au is not None and r_s_au is not None and temp_for_rxn is not None and temp_for_rxn > 0:
                    dS_rxn = (p_S_au - r_s_au) * J_PER_MOL_K_FROM_HARTREE_PER_K * CAL_PER_MOL_K_FROM_J_PER_MOL_K
                    TdS_rxn = ((p_S_au - r_s_au) * temp_for_rxn) * HARTREE_TO_KCAL_PER_MOL
                    # dH can be derived from dG and TdS
                    dH_rxn = dG_rxn + TdS_rxn; notes_rxn += f", ΔH from G,S"
                # Fallback for dH if S terms are missing
                elif p_H is not None and r_H_for_rxn is not None:
                     dH_rxn = (p_H - r_H_for_rxn) * HARTREE_TO_KCAL_PER_MOL; notes_rxn += f", ΔH({pH_type}-{rH_type})"
                else: dH_rxn = "N/A" # Cannot determine dH
            elif p_H is not None and r_H_for_rxn is not None: # Only H available
                 dH_rxn = (p_H - r_H_for_rxn) * HARTREE_TO_KCAL_PER_MOL; notes_rxn = f"ΔH(H_P-H_R)"
            
            insert_ch_results(f"R → P ({p_name_display})", dG_rxn, dH_rxn, -TdS_rxn if isinstance(TdS_rxn,float) else TdS_rxn, dS_rxn, "", notes_rxn)
            self.ch_results_tree.insert('', tk.END, values=("---","---","---","---","---","---","---"), tags=('separator',))

        ts_analysis_data = []
        for ts_detail in ts_details_for_calc:
            ts_path, ts_name_short, ts_data = ts_detail['path'], ts_detail['name'], ts_detail['data']
            if ts_data is None: # Skip if TS data failed to parse
                insert_ch_results(f"R → TS ({ts_name_short})", "ParseErr", "ParseErr", "ParseErr", "ParseErr", "", "TS data parsing failed")
                continue

            ts_G, tsG_type, ts_temp_eff = self._get_relevant_energy(ts_data, ['G'])
            ts_H, tsH_type, _ = self._get_relevant_energy(ts_data, ['H','E_el+ZPE','E_el'])
            ts_S_au, _, _ = self._get_s_and_ts(ts_data)
            r_G_for_act, rG_type, _ = self._get_relevant_energy(r_data, ['G']) # Use G of reactant for dG_act
            r_H_for_act, rH_type, _ = self._get_relevant_energy(r_data, ['H','E_el+ZPE','E_el']) # Use H of reactant for dH_act
            
            dG_act, dH_act, dS_act, TdS_act, notes_act = "N/A", "N/A", "N/A", "N/A", "Missing data"
            temp_for_act = ts_temp_eff if ts_temp_eff is not None else r_temp # Prioritize TS temp

            if ts_G is not None and r_G_for_act is not None:
                dG_act = (ts_G - r_G_for_act) * HARTREE_TO_KCAL_PER_MOL; notes_act = f"ΔG(G_TS-G_R)"
                if ts_S_au is not None and r_s_au is not None and temp_for_act is not None and temp_for_act > 0:
                    dS_act = (ts_S_au - r_s_au) * J_PER_MOL_K_FROM_HARTREE_PER_K * CAL_PER_MOL_K_FROM_J_PER_MOL_K
                    TdS_act = ((ts_S_au - r_s_au) * temp_for_act) * HARTREE_TO_KCAL_PER_MOL
                    dH_act = dG_act + TdS_act; notes_act += f", ΔH from G,S"
                elif ts_H is not None and r_H_for_act is not None:
                    dH_act = (ts_H - r_H_for_act) * HARTREE_TO_KCAL_PER_MOL; notes_act += f", ΔH({tsH_type}-{rH_type})"
                else: dH_act = "N/A"
            elif ts_H is not None and r_H_for_act is not None:
                dH_act = (ts_H - r_H_for_act) * HARTREE_TO_KCAL_PER_MOL; notes_act = f"ΔH(H_TS-H_R)"

            insert_ch_results(f"R → TS ({ts_name_short})", dG_act, dH_act, -TdS_act if isinstance(TdS_act,float) else TdS_act, dS_act, "", notes_act)
            if isinstance(dG_act, float): ts_analysis_data.append({'name': ts_name_short, 'dG_kcal': dG_act})

        if len(ts_analysis_data) >= 1 and r_temp is not None: # Ensure r_temp is valid
            rt_kcal = R_KCAL_MOL_K * r_temp
            if rt_kcal < 1e-9: messagebox.showwarning("CH Warning", "Reactant temp invalid or zero for ratios."); return
            
            min_dG_kcal_among_ts = min(ts['dG_kcal'] for ts in ts_analysis_data) if ts_analysis_data else 0
            boltzmann_factors = []; total_boltzmann_factor = 0
            for ts_res in ts_analysis_data:
                delta_delta_g = ts_res['dG_kcal'] - min_dG_kcal_among_ts; factor = 0.0
                try: factor = math.exp(-delta_delta_g / rt_kcal)
                except OverflowError: factor = 0.0 if delta_delta_g > 0 else float('inf') # Avoid large positive exp
                
                # Handle cases where factor might become inf due to identical dG to min_dG
                if math.isinf(factor) and delta_delta_g == 0: factor = 1.0 # Should be exp(0) = 1
                elif math.isinf(factor) and total_boltzmann_factor > 0 : factor = total_boltzmann_factor * 1e6 # Heuristic for large differences
                elif math.isinf(factor) : factor = 1e12 # Default large factor if it's the only one or first

                boltzmann_factors.append({'name': ts_res['name'], 'factor': factor})
                if not math.isinf(factor): total_boltzmann_factor += factor
            
            if math.isinf(total_boltzmann_factor) and len(boltzmann_factors)>0:
                num_infs = sum(1 for bf in boltzmann_factors if math.isinf(bf['factor']) or bf['factor'] > 1e10) # Count very large factors as inf for ratio
                for bf_item in boltzmann_factors: 
                    bf_item['percentage'] = 100.0 / num_infs if num_infs > 0 and (math.isinf(bf_item['factor']) or bf_item['factor'] > 1e10) else 0.0
            elif total_boltzmann_factor > 1e-9:
                for bf_item in boltzmann_factors: bf_item['percentage'] = (bf_item['factor'] / total_boltzmann_factor) * 100
            else: # All factors are essentially zero or only one TS
                 for i, bf_item in enumerate(boltzmann_factors): bf_item['percentage'] = 100.0 if len(boltzmann_factors) == 1 and i==0 else 0.0
            
            all_ch_items = self.ch_results_tree.get_children('')
            for item_id in all_ch_items:
                current_values = list(self.ch_results_tree.item(item_id, 'values'))
                pathway_name = current_values[0]
                for bf_item in boltzmann_factors:
                    if bf_item['name'] in pathway_name and "R → TS" in pathway_name:
                        current_values[5] = f"{bf_item['percentage']:.1f}%"; self.ch_results_tree.item(item_id, values=tuple(current_values)); break
        self.apply_treeview_row_tags(self.ch_results_tree)

    def plot_curtin_hammett_profile(self): 
        if not MATPLOTLIB_AVAILABLE:
            messagebox.showerror("Plot Error", "Matplotlib not installed.")
            return

        r_name_display = self.ch_reactant_var.get()
        p_name_display = self.ch_product_var.get()
        ch_ts_display_names = list(self.ch_ts_listbox.get(0, tk.END))

        r_path = self._get_ch_filepath_from_display_name(r_name_display)
        if not r_path or r_path not in self.loaded_files_data or self.loaded_files_data[r_path] is None:
            messagebox.showerror("CH Plot Error", "Reactant data missing or invalid.")
            return
        
        r_data = self.loaded_files_data[r_path]
        r_energy_eh, r_energy_type, _ = self._get_relevant_energy(r_data, ['G', 'H', 'E_el+ZPE', 'E_el']) 
        if r_energy_eh is None:
            messagebox.showerror("CH Plot Error", f"Could not retrieve a suitable energy for reactant: {r_name_display}")
            return

        plot_points = []
        plot_labels = []
        plot_energy_types = [] 

        plot_points.append(r_energy_eh)
        plot_labels.append(f"R ({os.path.basename(r_path).split('.')[0]})")
        plot_energy_types.append(r_energy_type)

        for ts_name_disp in ch_ts_display_names:
            ts_path = self._get_ch_filepath_from_display_name(ts_name_disp)
            if ts_path and ts_path in self.loaded_files_data and self.loaded_files_data[ts_path]:
                ts_data = self.loaded_files_data[ts_path]
                ts_energy_eh, ts_energy_type, _ = self._get_relevant_energy(ts_data, ['G', 'H', 'E_el+ZPE', 'E_el'])
                if ts_energy_eh is not None:
                    plot_points.append(ts_energy_eh)
                    plot_labels.append(f"TS ({os.path.basename(ts_path).split('.')[0]})")
                    plot_energy_types.append(ts_energy_type)
                else:
                    messagebox.showwarning("CH Plot Warning", f"Could not retrieve energy for TS: {ts_name_disp}. Skipping in plot.")
            else:
                 messagebox.showwarning("CH Plot Warning", f"Data for TS: {ts_name_disp} not found. Skipping in plot.")

        p_energy_eh = None
        if p_name_display and p_name_display != "(None)":
            p_path = self._get_ch_filepath_from_display_name(p_name_display)
            if p_path and p_path in self.loaded_files_data and self.loaded_files_data[p_path]:
                p_data = self.loaded_files_data[p_path]
                p_energy_eh, p_energy_type, _ = self._get_relevant_energy(p_data, ['G', 'H', 'E_el+ZPE', 'E_el'])
                if p_energy_eh is not None:
                    plot_points.append(p_energy_eh)
                    plot_labels.append(f"P ({os.path.basename(p_path).split('.')[0]})")
                    plot_energy_types.append(p_energy_type)
                else:
                    messagebox.showwarning("CH Plot Warning", f"Could not retrieve energy for product: {p_name_display}. Skipping in plot.")
            else:
                messagebox.showwarning("CH Plot Warning", f"Data for product: {p_name_display} not found. Skipping in plot.")
        
        if len(plot_points) < 2 : 
            messagebox.showinfo("CH Plot Info", "Not enough data points (reactant and at least one TS) to plot a profile.")
            return

        relative_energies_kcal = [(e - r_energy_eh) * HARTREE_TO_KCAL_PER_MOL for e in plot_points]
        
        y_label = "Gibbs Free Energy (G)" 
        if r_energy_type != 'G':
            y_label = f"{r_energy_type.replace('_el', 'E_el')} (used for Reactant)"

        line_style = self.ch_plot_line_style_var.get()
        line_width_str = self.ch_plot_line_width_var.get()
        line_color = self.ch_plot_line_color_var.get()

        self._plot_profile_matplotlib(
            x_coords=range(len(plot_labels)),
            y_values_kcal=relative_energies_kcal,
            x_labels=plot_labels,
            y_label_main=y_label, 
            title_text="Curtin-Hammett Energy Profile",
            ref_label_for_plot=f"Reactant ({os.path.basename(r_path).split('.')[0]})",
            types_plotted_for_annotations=plot_energy_types,
            enable_break=self.enable_ch_axis_break_var.get(),
            break_from_str=self.ch_axis_break_from_var.get(),
            break_to_str=self.ch_axis_break_to_var.get(),
            line_style=line_style,
            line_width_str=line_width_str,
            line_color=line_color,
            is_allxyz_plot=False 
        )


    # --- Methods for TAB 2 (TS Estimation) ---
    def move_file_in_list(self, direction): 
        selected_indices = self.file_listbox.curselection()
        if not selected_indices:
            messagebox.showinfo("Reorder Info", "Please select a file in the list to move.")
            return
        
        current_index = selected_indices[0]
        
        if direction == "up":
            if current_index == 0: return
            new_index = current_index - 1
        elif direction == "down":
            if current_index == self.file_listbox.size() - 1: return
            new_index = current_index + 1
        else:
            return

        path_to_move = self.loaded_file_paths.pop(current_index)
        self.loaded_file_paths.insert(new_index, path_to_move)
        
        # Re-populate the listbox based on the new order of loaded_file_paths
        self.file_listbox.delete(0, tk.END)
        for f_path in self.loaded_file_paths:
            filename = os.path.basename(f_path)
            self.file_listbox.insert(tk.END, filename)
            # Re-apply red color if parsing failed for this file
            if f_path not in self.loaded_files_data or self.loaded_files_data[f_path] is None:
                try: 
                    idx_in_listbox = self.file_listbox.get(0, tk.END).index(filename)
                    self.file_listbox.itemconfig(idx_in_listbox, {'fg': 'red'})
                except ValueError:
                    pass 

        self.file_listbox.selection_clear(0, tk.END) 
        self.file_listbox.selection_set(new_index)
        self.file_listbox.activate(new_index)
        self.file_listbox.see(new_index)
        
        self.on_file_select_from_listbox() 

        self.update_plot_reference_selector()
        self.update_ch_selectors()
        self.check_plot_button_state()


    def load_allxyz_files(self):
        filepaths_tuple = filedialog.askopenfilenames(title="Select .allxyz or multi-frame .xyz File(s)",filetypes=(("XYZ trajectory files", "*.allxyz *.xyz"), ("All files", "*.*")))
        if not filepaths_tuple: return
        self.allxyz_file_paths = list(filepaths_tuple); self.allxyz_files_data.clear(); self.allxyz_listbox.delete(0, tk.END)
        self.current_selected_allxyz_filepath = None; self.allxyz_filepath_label.config(text="No .allxyz file selected")
        for item in self.allxyz_details_tree.get_children(): self.allxyz_details_tree.delete(item)
        self.highest_e_info_label.config(text="Highest Energy Frame: N/A"); self.save_highest_e_button.config(state=tk.DISABLED); self.plot_allxyz_button.config(state=tk.DISABLED)
        self.highest_energy_frame_data = None
        for filepath in self.allxyz_file_paths:
            filename = os.path.basename(filepath); self.allxyz_listbox.insert(tk.END, filename)
            frames = parse_multiframe_xyz(filepath) # This now has encoding fallback
            if frames: self.allxyz_files_data[filepath] = frames
            else: 
                try:
                    idx_to_color = self.allxyz_listbox.get(0, tk.END).index(filename)
                    self.allxyz_listbox.itemconfig(idx_to_color, {'fg': 'red'})
                except ValueError:
                    pass
        if self.allxyz_listbox.size() > 0: self.allxyz_listbox.selection_set(0); self.allxyz_listbox.event_generate("<<ListboxSelect>>")

    def on_allxyz_file_select(self, event=None):
        widget = event.widget if event else self.allxyz_listbox; selected_indices = widget.curselection()
        if not selected_indices: return
        try: self.current_selected_allxyz_filepath = self.allxyz_file_paths[selected_indices[0]]
        except IndexError: self.allxyz_filepath_label.config(text="Error selecting .allxyz file"); return
        selected_filename = os.path.basename(self.current_selected_allxyz_filepath); self.allxyz_filepath_label.config(text=f"Selected: {selected_filename}")
        for item in self.allxyz_details_tree.get_children(): self.allxyz_details_tree.delete(item)
        self.highest_e_info_label.config(text="Highest Energy Frame: N/A"); self.save_highest_e_button.config(state=tk.DISABLED)
        self.plot_allxyz_button.config(state=tk.DISABLED); self.highest_energy_frame_data = None
        
        if self.current_selected_allxyz_filepath in self.allxyz_files_data and self.allxyz_files_data[self.current_selected_allxyz_filepath] is not None:
            frames = self.allxyz_files_data[self.current_selected_allxyz_filepath]
            if not frames: self.allxyz_details_tree.insert('', tk.END, values=("Info", "No valid frames.", "")); return
            
            # Filter out frames where energy_eh is None for min/max operations
            valid_energy_frames = [f for f in frames if f['energy_eh'] is not None]
            if not valid_energy_frames:
                # Populate with N/A if no frames have energy
                for frame_data in frames:
                     self.allxyz_details_tree.insert('', tk.END, values=(frame_data['frame_index'] + 1, "N/A", "N/A"))
                self.highest_e_info_label.config(text="Highest Energy Frame: N/A (No energy data)")
                return # No valid energies to plot or find highest

            min_e = min(f['energy_eh'] for f in valid_energy_frames)
            highest_e = -float('inf'); temp_highest_frame_data = None
            
            for frame_data in frames: # Iterate all frames for display
                energy_eh = frame_data['energy_eh']; rel_e_kcal_str = "N/A"
                tags_to_apply = ()
                if energy_eh is not None: # Process only if energy exists for this frame
                    rel_e_kcal_str = f"{(energy_eh - min_e) * HARTREE_TO_KCAL_PER_MOL:.2f}"
                    if energy_eh > highest_e: 
                        highest_e = energy_eh; temp_highest_frame_data = frame_data
                    if temp_highest_frame_data and frame_data['frame_index'] == temp_highest_frame_data['frame_index']:
                        tags_to_apply = ('highest_e',)
                
                self.allxyz_details_tree.insert('', tk.END, values=(frame_data['frame_index'] + 1, f"{energy_eh:.8f}" if energy_eh is not None else "N/A", rel_e_kcal_str), tags=tags_to_apply)
            
            if temp_highest_frame_data:
                self.highest_energy_frame_data = temp_highest_frame_data
                self.highest_e_info_label.config(text=f"Highest E: Frame {temp_highest_frame_data['frame_index']+1}, E = {temp_highest_frame_data['energy_eh']:.8f} Eh")
                self.save_highest_e_button.config(state=tk.NORMAL)
            
            if MATPLOTLIB_AVAILABLE and valid_energy_frames: self.plot_allxyz_button.config(state=tk.NORMAL)
        else: self.allxyz_details_tree.insert('', tk.END, values=("Error", "Data not loaded or parsing failed.", ""))


    def save_highest_energy_geometry(self):
        if not self.highest_energy_frame_data: messagebox.showerror("Error", "No highest E frame."); return
        default_filename = f"ts_guess_from_{os.path.basename(self.highest_energy_frame_data['filepath'])}_frame{self.highest_energy_frame_data['frame_index']+1}.xyz"
        save_path = filedialog.asksaveasfilename(title="Save Highest Energy Geometry as XYZ", defaultextension=".xyz", initialfile=default_filename, filetypes=(("XYZ files", "*.xyz"), ("All files", "*.*")))
        if not save_path: return
        try:
            with open(save_path, 'w', encoding='utf-8') as f: # XYZ usually safe with utf-8
                f.write(f"{self.highest_energy_frame_data['num_atoms']}\n")
                f.write(f"{self.highest_energy_frame_data['comment_line']} (Frame {self.highest_energy_frame_data['frame_index']+1})\n")
                for coord_line in self.highest_energy_frame_data['coordinates']: f.write(coord_line + "\n")
            messagebox.showinfo("Success", f"Saved to:\n{save_path}")
        except Exception as e: messagebox.showerror("Save Error", f"Could not save: {e}")

    def plot_allxyz_profile(self):
        if not MATPLOTLIB_AVAILABLE: messagebox.showerror("Plot Error", "Matplotlib not installed."); return
        if not self.current_selected_allxyz_filepath or self.current_selected_allxyz_filepath not in self.allxyz_files_data or self.allxyz_files_data[self.current_selected_allxyz_filepath] is None: 
            messagebox.showerror("Plot Error", "No .allxyz file selected or data not loaded/parsed."); return
        
        frames = self.allxyz_files_data[self.current_selected_allxyz_filepath]; 
        valid_frames = [f for f in frames if f['energy_eh'] is not None]
        if not valid_frames: messagebox.showinfo("Plot Info", "No energy data in selected .allxyz file to plot."); return
        
        frame_indices = [f['frame_index'] + 1 for f in valid_frames]; energies_eh = [f['energy_eh'] for f in valid_frames]
        min_e = min(energies_eh); rel_energies_kcal = [(e - min_e) * HARTREE_TO_KCAL_PER_MOL for e in energies_eh]
        
        # Find the frame index (1-based) corresponding to the min_e for labeling
        min_e_frame_index_for_label = frame_indices[energies_eh.index(min_e)] if energies_eh else "N/A"

        line_style = self.allxyz_plot_line_style_var.get()
        line_width_str = self.allxyz_plot_line_width_var.get()
        line_color = self.allxyz_plot_line_color_var.get()

        self._plot_profile_matplotlib(
            x_coords=frame_indices,
            y_values_kcal=rel_energies_kcal,
            x_labels=[str(fi) for fi in frame_indices], # Use frame indices as x-labels
            y_label_main="Energy", # Generic y-label
            title_text=f"Energy Profile for {os.path.basename(self.current_selected_allxyz_filepath)}",
            ref_label_for_plot=f"Frame {min_e_frame_index_for_label}, E={min_e:.6f} Eh",
            types_plotted_for_annotations=None, # No specific types like G, H for allxyz frames
            enable_break=self.enable_allxyz_axis_break_var.get(),
            break_from_str=self.allxyz_axis_break_from_var.get(),
            break_to_str=self.allxyz_axis_break_to_var.get(),
            line_style=line_style,
            line_width_str=line_width_str,
            line_color=line_color,
            is_allxyz_plot=True
        )


    def copy_allxyz_table_to_clipboard(self):
        if not self.allxyz_details_tree.get_children() or \
           (len(self.allxyz_details_tree.get_children()) == 1 and self.allxyz_details_tree.item(self.allxyz_details_tree.get_children()[0], 'values')[0] in ["Error", "Info"]):
            messagebox.showinfo("Clipboard", ".allxyz table is empty or contains no data. Nothing to copy.")
            return
        header = ["Frame #", "Energy (Eh)", "Rel. E (kcal/mol)"]
        tsv_data = ["\t".join(header)]
        for item_id in self.allxyz_details_tree.get_children():
            values = self.allxyz_details_tree.item(item_id, 'values')
            row_values = [str(v) if v is not None else "" for v in values]
            tsv_data.append("\t".join(row_values))
        full_tsv_string = "\n".join(tsv_data)
        try:
            self.master.clipboard_clear()
            self.master.clipboard_append(full_tsv_string)
            messagebox.showinfo("Clipboard", ".allxyz frame table data copied to clipboard as TSV.")
        except tk.TclError:
            messagebox.showerror("Clipboard Error", "Could not access the clipboard.")

    # --- Consolidated Table Methods ---
    def open_consolidated_table_window(self):
        if not self.loaded_file_paths:
            messagebox.showinfo("Consolidated Table", "No files loaded to display in the table.")
            return

        if self.consolidated_table_window and self.consolidated_table_window.winfo_exists():
            self.consolidated_table_window.lift()
            return

        self.consolidated_table_window = tk.Toplevel(self.master)
        self.consolidated_table_window.title("Consolidated Energy Table")
        self.consolidated_table_window.geometry("700x500")

        options_frame = ttk.LabelFrame(self.consolidated_table_window, text="Table Options", padding=10)
        options_frame.pack(padx=10, pady=10, fill=tk.X)

        ttk.Label(options_frame, text="Energy Parameter:").grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.consolidated_param_combo = ttk.Combobox(options_frame, textvariable=self.consolidated_param_var,
                                                     values=["Gibbs Free Energy (G)", "Enthalpy (H)", "E_el + ZPVE", "Electronic Energy (E_el)"],
                                                     state="readonly", width=25)
        self.consolidated_param_combo.grid(row=0, column=1, padx=5, pady=5, sticky=tk.EW)
        self.consolidated_param_var.set("Gibbs Free Energy (G)")
        self.consolidated_param_combo.bind("<<ComboboxSelected>>", lambda e: self.populate_consolidated_table_view())
        
        ttk.Label(options_frame, text="Display Unit:").grid(row=0, column=2, padx=5, pady=5, sticky=tk.W)
        self.consolidated_unit_combo = ttk.Combobox(options_frame, textvariable=self.consolidated_unit_var,
                                                    values=["kcal/mol", "kJ/mol", "Eh", "eV"],
                                                    state="readonly", width=15)
        self.consolidated_unit_combo.grid(row=0, column=3, padx=5, pady=5, sticky=tk.EW)
        self.consolidated_unit_var.set("kcal/mol")
        self.consolidated_unit_combo.bind("<<ComboboxSelected>>", lambda e: self.populate_consolidated_table_view())
        
        options_frame.grid_columnconfigure(1, weight=1)
        options_frame.grid_columnconfigure(3, weight=1)

        tree_frame = ttk.Frame(self.consolidated_table_window)
        tree_frame.pack(padx=10, pady=(0,10), fill=tk.BOTH, expand=True)

        self.consolidated_tree = ttk.Treeview(tree_frame, columns=("file", "value", "unit"), show="headings", selectmode="extended")
        self.consolidated_tree.heading("file", text="File", anchor=tk.W)
        self.consolidated_tree.heading("value", text="Energy Value", anchor=tk.E)
        self.consolidated_tree.heading("unit", text="Unit", anchor=tk.W)
        self.consolidated_tree.column("file", width=300, stretch=tk.YES, anchor=tk.W)
        self.consolidated_tree.column("value", width=150, stretch=tk.NO, anchor=tk.E)
        self.consolidated_tree.column("unit", width=100, stretch=tk.NO, anchor=tk.W)
        
        vsb = ttk.Scrollbar(tree_frame, orient="vertical", command=self.consolidated_tree.yview)
        hsb = ttk.Scrollbar(tree_frame, orient="horizontal", command=self.consolidated_tree.xview)
        self.consolidated_tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        vsb.pack(side='right', fill='y'); hsb.pack(side='bottom', fill='x')
        self.consolidated_tree.pack(side='left', fill='both', expand=True)
        self.apply_treeview_row_tags(self.consolidated_tree) # Apply odd/even row styling

        button_frame = ttk.Frame(self.consolidated_table_window)
        button_frame.pack(pady=5, fill=tk.X, padx=10)
        ttk.Button(button_frame, text="Copy Table (TSV)", command=self.copy_consolidated_table_to_clipboard).pack(side=tk.LEFT)
        ttk.Button(button_frame, text="Refresh Table", command=self.populate_consolidated_table_view).pack(side=tk.LEFT, padx=5)


        self.populate_consolidated_table_view()
        self.consolidated_table_window.protocol("WM_DELETE_WINDOW", self._on_consolidated_close)

    def _on_consolidated_close(self):
        if self.consolidated_table_window:
            self.consolidated_table_window.destroy()
            self.consolidated_table_window = None


    def populate_consolidated_table_view(self):
        if not self.consolidated_table_window or not self.consolidated_table_window.winfo_exists():
            return # Window doesn't exist or has been closed

        for item in self.consolidated_tree.get_children():
            self.consolidated_tree.delete(item)

        selected_param_str = self.consolidated_param_var.get()
        selected_unit_str = self.consolidated_unit_var.get()

        param_priority_map = {
            "Gibbs Free Energy (G)": ['G', 'H', 'E_el+ZPE', 'E_el'],
            "Enthalpy (H)": ['H', 'E_el+ZPE', 'E_el'],
            "E_el + ZPVE": ['E_el+ZPE', 'E_el'],
            "Electronic Energy (E_el)": ['E_el']
        }
        priority = param_priority_map.get(selected_param_str, ['G', 'H', 'E_el+ZPE', 'E_el'])


        for path in self.loaded_file_paths:
            data = self.loaded_files_data.get(path)
            filename = os.path.basename(path)
            energy_eh, energy_type_found, _ = self._get_relevant_energy(data, priority)
            
            display_value = "N/A"
            actual_unit_str = ""

            if energy_eh is not None:
                actual_unit_str = selected_unit_str # Store the target unit
                if selected_unit_str == "kcal/mol":
                    converted_value = energy_eh * HARTREE_TO_KCAL_PER_MOL
                    display_value = f"{converted_value:.2f}"
                elif selected_unit_str == "kJ/mol":
                    converted_value = energy_eh * HARTREE_TO_KJ_PER_MOL
                    display_value = f"{converted_value:.2f}"
                elif selected_unit_str == "eV":
                    converted_value = energy_eh * HARTREE_TO_EV
                    display_value = f"{converted_value:.4f}"
                elif selected_unit_str == "Eh":
                    display_value = f"{energy_eh:.8f}"
                else: # Should not happen with combobox
                    display_value = f"{energy_eh:.8f}" 
                    actual_unit_str = "Eh (default)"
            else: # energy_eh is None
                if data is None: # Parsing failed completely
                    display_value = "Parse Err"
                else: # Parsing succeeded but specific energy not found
                    display_value = f"N/A ({energy_type_found})" # energy_type_found might be N/A if nothing was found
                actual_unit_str = "-"


            self.consolidated_tree.insert('', tk.END, values=(filename, display_value, actual_unit_str))
        self.apply_treeview_row_tags(self.consolidated_tree)


    def copy_consolidated_table_to_clipboard(self):
        if not self.consolidated_table_window or not self.consolidated_table_window.winfo_exists() or \
           not self.consolidated_tree.get_children():
            messagebox.showinfo("Clipboard", "Consolidated table is empty or not available.")
            return
        
        header = ["File", "Energy Value", "Unit"]
        tsv_data = ["\t".join(header)]
        for item_id in self.consolidated_tree.get_children():
            values = self.consolidated_tree.item(item_id, 'values')
            tsv_data.append("\t".join([str(v) if v is not None else "" for v in values]))
        
        full_tsv_string = "\n".join(tsv_data)
        try:
            self.master.clipboard_clear()
            self.master.clipboard_append(full_tsv_string)
            messagebox.showinfo("Clipboard", "Consolidated table data copied to clipboard (TSV).", parent=self.consolidated_table_window)
        except tk.TclError:
            messagebox.showerror("Clipboard Error", "Could not access the clipboard.", parent=self.consolidated_table_window)

    # --- Step Delta Calculation Methods ---
    def open_step_deltas_window(self):
        if not self.loaded_file_paths:
            messagebox.showinfo("Step Δ Values", "No files loaded to define steps.")
            return

        if self.step_delta_window and self.step_delta_window.winfo_exists():
            self.step_delta_window.lift()
            return

        self.step_delta_window = tk.Toplevel(self.master)
        self.step_delta_window.title("Stepwise Thermodynamic Changes")
        self.step_delta_window.geometry("850x500")

        controls_frame = ttk.LabelFrame(self.step_delta_window, text="Define Step", padding=10)
        controls_frame.pack(padx=10, pady=10, fill=tk.X)

        ttk.Label(controls_frame, text="Reactant for Step:").grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
        self.step_reactant_combo = ttk.Combobox(controls_frame, textvariable=self.step_reactant_var, state="readonly", width=35)
        self.step_reactant_combo.grid(row=0, column=1, padx=5, pady=5, sticky=tk.EW)

        ttk.Label(controls_frame, text="Product for Step:").grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)
        self.step_product_combo = ttk.Combobox(controls_frame, textvariable=self.step_product_var, state="readonly", width=35)
        self.step_product_combo.grid(row=1, column=1, padx=5, pady=5, sticky=tk.EW)
        
        controls_frame.grid_columnconfigure(1, weight=1)

        add_button = ttk.Button(controls_frame, text="Add Step to Table", command=self._calculate_and_add_step_to_table)
        add_button.grid(row=0, column=2, rowspan=2, padx=10, pady=5, sticky=tk.NS)
        
        results_frame = ttk.LabelFrame(self.step_delta_window, text="Calculated Step Values", padding=10)
        results_frame.pack(padx=10, pady=5, fill=tk.BOTH, expand=True)

        cols = ("step", "delta_g", "delta_h", "delta_s", "neg_tds", "temp_r", "temp_p")
        self.step_delta_tree = ttk.Treeview(results_frame, columns=cols, show="headings")
        self.step_delta_tree.heading("step", text="Step (R → P)", anchor=tk.W)
        self.step_delta_tree.heading("delta_g", text="ΔG (kcal/mol)", anchor=tk.E)
        self.step_delta_tree.heading("delta_h", text="ΔH (kcal/mol)", anchor=tk.E)
        self.step_delta_tree.heading("delta_s", text="ΔS (cal/mol·K)", anchor=tk.E)
        self.step_delta_tree.heading("neg_tds", text="-TΔS (kcal/mol)", anchor=tk.E)
        self.step_delta_tree.heading("temp_r", text="T_R (K)", anchor=tk.E)
        self.step_delta_tree.heading("temp_p", text="T_P (K)", anchor=tk.E)

        self.step_delta_tree.column("step", width=250, stretch=tk.YES, anchor=tk.W)
        for col in cols[1:]: self.step_delta_tree.column(col, width=100, stretch=tk.NO, anchor=tk.E)
        
        vsb_step = ttk.Scrollbar(results_frame, orient="vertical", command=self.step_delta_tree.yview)
        hsb_step = ttk.Scrollbar(results_frame, orient="horizontal", command=self.step_delta_tree.xview)
        self.step_delta_tree.configure(yscrollcommand=vsb_step.set, xscrollcommand=hsb_step.set)
        vsb_step.pack(side='right', fill='y'); hsb_step.pack(side='bottom', fill='x')
        self.step_delta_tree.pack(fill=tk.BOTH, expand=True)
        self.apply_treeview_row_tags(self.step_delta_tree)

        bottom_buttons_frame = ttk.Frame(self.step_delta_window)
        bottom_buttons_frame.pack(pady=5, fill=tk.X, padx=10)
        ttk.Button(bottom_buttons_frame, text="Copy Table (TSV)", command=self._copy_step_deltas_table).pack(side=tk.LEFT)
        ttk.Button(bottom_buttons_frame, text="Clear Table", command=self._clear_step_deltas_table).pack(side=tk.LEFT, padx=5)


        self._populate_step_deltas_combos()
        self._display_step_deltas_in_treeview() # Display any previously stored data
        self.step_delta_window.protocol("WM_DELETE_WINDOW", self._on_step_deltas_close)

    def _populate_step_deltas_combos(self):
        if not self.step_delta_window or not self.step_delta_window.winfo_exists(): return
        
        valid_files = [os.path.basename(p) for p in self.loaded_file_paths if p in self.loaded_files_data and self.loaded_files_data[p]]
        valid_files.sort()
        
        self.step_reactant_combo['values'] = valid_files
        self.step_product_combo['values'] = valid_files
        if valid_files:
            if not self.step_reactant_var.get() or self.step_reactant_var.get() not in valid_files:
                self.step_reactant_var.set(valid_files[0])
            if not self.step_product_var.get() or self.step_product_var.get() not in valid_files:
                 self.step_product_var.set(valid_files[0] if len(valid_files) == 1 else (valid_files[1] if len(valid_files) > 1 else ""))
        else:
            self.step_reactant_var.set("")
            self.step_product_var.set("")


    def _calculate_and_add_step_to_table(self):
        r_name_disp = self.step_reactant_var.get()
        p_name_disp = self.step_product_var.get()

        if not r_name_disp or not p_name_disp:
            messagebox.showerror("Input Error", "Please select both a reactant and a product for the step.", parent=self.step_delta_window)
            return

        r_path = self._get_ch_filepath_from_display_name(r_name_disp) # Reusing this helper
        p_path = self._get_ch_filepath_from_display_name(p_name_disp)

        if not r_path or r_path not in self.loaded_files_data or self.loaded_files_data[r_path] is None:
            messagebox.showerror("Data Error", f"Data for reactant '{r_name_disp}' not found or invalid.", parent=self.step_delta_window)
            return
        if not p_path or p_path not in self.loaded_files_data or self.loaded_files_data[p_path] is None:
            messagebox.showerror("Data Error", f"Data for product '{p_name_disp}' not found or invalid.", parent=self.step_delta_window)
            return
        
        if r_path == p_path:
            messagebox.showwarning("Input Error", "Reactant and Product for the step cannot be the same file.", parent=self.step_delta_window)
            return

        r_data = self.loaded_files_data[r_path]
        p_data = self.loaded_files_data[p_path]

        # Get G, H, S, T for Reactant
        g_r, _, t_r = self._get_relevant_energy(r_data, ['G'])
        h_r, _, _ = self._get_relevant_energy(r_data, ['H']) # Temp from G is fine
        s_r_au, ts_r_eh, t_r_s = self._get_s_and_ts(r_data) # t_r_s should match t_r if thermo block consistent

        # Get G, H, S, T for Product
        g_p, _, t_p = self._get_relevant_energy(p_data, ['G'])
        h_p, _, _ = self._get_relevant_energy(p_data, ['H'])
        s_p_au, ts_p_eh, t_p_s = self._get_s_and_ts(p_data)

        # Default to 298.15 K if temperatures are None
        t_r_eff = t_r if t_r is not None else 298.15
        t_p_eff = t_p if t_p is not None else 298.15
        
        # Calculations
        delta_g_kcal, delta_h_kcal, delta_s_cal, neg_tds_kcal = "N/A", "N/A", "N/A", "N/A"

        if g_r is not None and g_p is not None:
            delta_g_kcal = (g_p - g_r) * HARTREE_TO_KCAL_PER_MOL
        
        if h_r is not None and h_p is not None:
            delta_h_kcal = (h_p - h_r) * HARTREE_TO_KCAL_PER_MOL

        if s_r_au is not None and s_p_au is not None:
            delta_s_au = s_p_au - s_r_au
            delta_s_cal = delta_s_au * J_PER_MOL_K_FROM_HARTREE_PER_K * CAL_PER_MOL_K_FROM_J_PER_MOL_K
            # Use reactant's temperature for TdS term, or product's if reactant's is missing
            temp_for_tds = t_r_eff 
            if t_r is None and t_p is not None: temp_for_tds = t_p_eff

            neg_tds_kcal = - (temp_for_tds * delta_s_cal / 1000.0) # Convert cal to kcal
        
        # Store and display
        step_label = f"{r_name_disp.split('.')[0]} → {p_name_disp.split('.')[0]}"
        self.step_delta_results_data.append((
            step_label, delta_g_kcal, delta_h_kcal, delta_s_cal, neg_tds_kcal, 
            f"{t_r_eff:.2f}" if t_r is not None else "N/A", 
            f"{t_p_eff:.2f}" if t_p is not None else "N/A"
        ))
        self._display_step_deltas_in_treeview()

    def _display_step_deltas_in_treeview(self):
        if not self.step_delta_window or not self.step_delta_window.winfo_exists(): return
        for item in self.step_delta_tree.get_children(): self.step_delta_tree.delete(item)
        
        for row_data in self.step_delta_results_data:
            formatted_row = []
            for i, val in enumerate(row_data):
                if isinstance(val, float):
                    if i in [1,2,4]: # dG, dH, -TdS (kcal/mol)
                        formatted_row.append(f"{val:.2f}")
                    elif i == 3: # dS (cal/molK)
                         formatted_row.append(f"{val:.2f}")
                    else: # Should not happen for float, but as fallback
                        formatted_row.append(str(val))

                else: # String or N/A
                    formatted_row.append(str(val))
            self.step_delta_tree.insert('', tk.END, values=tuple(formatted_row))
        self.apply_treeview_row_tags(self.step_delta_tree)

    def _clear_step_deltas_table(self):
        self.step_delta_results_data.clear()
        self._display_step_deltas_in_treeview()


    def _copy_step_deltas_table(self):
        if not self.step_delta_window or not self.step_delta_window.winfo_exists() or \
           not self.step_delta_tree.get_children():
            messagebox.showinfo("Clipboard", "Step Deltas table is empty.", parent=self.step_delta_window)
            return
        
        header = ["Step (R → P)", "ΔG (kcal/mol)", "ΔH (kcal/mol)", "ΔS (cal/mol·K)", "-TΔS (kcal/mol)", "T_R (K)", "T_P (K)"]
        tsv_data = ["\t".join(header)]
        for item_id in self.step_delta_tree.get_children():
            values = self.step_delta_tree.item(item_id, 'values')
            tsv_data.append("\t".join([str(v) if v is not None else "" for v in values]))
        
        full_tsv_string = "\n".join(tsv_data)
        try:
            self.master.clipboard_clear()
            self.master.clipboard_append(full_tsv_string)
            messagebox.showinfo("Clipboard", "Step Deltas table data copied to clipboard (TSV).", parent=self.step_delta_window)
        except tk.TclError:
            messagebox.showerror("Clipboard Error", "Could not access the clipboard.", parent=self.step_delta_window)

    def _on_step_deltas_close(self):
        if self.step_delta_window:
            self.step_delta_window.destroy()
            self.step_delta_window = None


if __name__ == "__main__":
    root = tk.Tk()
    app = ThermoApp(root)
    root.mainloop()
