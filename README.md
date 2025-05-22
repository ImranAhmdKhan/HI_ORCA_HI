# HI_ORCA_HI
# ORCA Thermo & TS Tools

## Overview

The "ORCA Thermo & TS Tools" is a Python-based graphical user interface (GUI) application designed to assist computational chemists in parsing, analyzing, and visualizing data from ORCA output files (`.out`, `.log`). It also includes features for analyzing multi-frame XYZ trajectory files (often named `.allxyz`) for transition state (TS) estimation and for performing Curtin-Hammett analysis.

The application is built using Tkinter for the GUI and can optionally use `spiffyplot` (if available on the system) or Matplotlib for generating energy profile diagrams.

## Features

* **File Loading & Management:**
    * Load multiple ORCA output files (`.out`, `.log`) and XYZ files (`.xyz`, `.allxyz`).
    * Files are listed, and their order can be rearranged for reaction profile plotting.
    * Files that fail to parse or are missing key data are highlighted.
* **Tab 1: Thermochemistry & Reaction Analysis**
    * **Detailed Thermochemistry Display:**
        * Shows detailed thermochemical parameters (Electronic Energy, ZPE, Thermal Corrections for U, H, G, T\*S terms, Total Entropy) for a selected file.
        * Allows selection of energy units for display (Eh, kcal/mol, kJ/mol, eV).
        * Users can toggle the visibility of different thermochemistry sections.
    * **Reaction Energy Profile Plotting:**
        * Generates energy profile diagrams for the sequence of loaded files.
        * Users can select the energy type to plot (G, H, E+ZPE, E\_el).
        * Allows selection of a reference file or uses the overall lowest energy structure as the zero-point for relative energies.
        * Customizable plot color.
        * Outputs plots as SVG files using `spiffyplot` (if available).
    * **Curtin-Hammett Analysis:**
        * Select a reactant, one or more transition states, and an optional product from the loaded files.
        * Calculates and displays $\Delta G^{\ddagger}$, $\Delta H^{\ddagger}$, $-T\Delta S^{\ddagger}$, and $\Delta S^{\ddagger}$ for each R $\rightarrow$ TS pathway.
        * Calculates $\Delta G_{rxn}$, $\Delta H_{rxn}$, etc., if a product is specified.
        * Determines product ratios based on the Boltzmann distribution over the transition state energies.
        * Option to plot the Curtin-Hammett energy profile (G vs. reaction coordinate) with customizable color.
    * **Stepwise Thermodynamic Changes ($\Delta G, \Delta H, \Delta S$):**
        * A dedicated module to calculate $\Delta G$, $\Delta H$, and $\Delta S$ between any two selected loaded files (Reactant $\rightarrow$ Product for a step).
        * Displays results in a table, including temperatures of the reactant and product species.
    * **Consolidated Energy Table:**
        * View a table of a specific energy parameter (G, H, E+ZPE, E\_el) for all loaded files in a chosen unit.
    * **Clipboard Functionality:** Copy displayed tables (thermo details, CH results, consolidated table, step deltas) to the clipboard in TSV format.
* **Tab 2: TS Estimation (.allxyz)**
    * Load multi-frame XYZ files (typically `.allxyz` from ORCA NEB or scan calculations).
    * Displays energy for each frame and relative energy (kcal/mol) with respect to the lowest energy frame.
    * Identifies and highlights the highest energy frame, often a good initial guess for a transition state.
    * **Save Highest Energy Geometry:** Save the coordinates of the highest energy frame as a new `.xyz` file.
    * **Plot Profile:** Generates an energy profile (Energy vs. Frame Number) for the loaded `.allxyz` file using `spiffyplot` with customizable color.
    * **Clipboard Functionality:** Copy the frame table to the clipboard in TSV format.
* **Tab 3: Orbital Analysis**
    * **HOMO/LUMO Display:** For the currently selected file (from Tab 1), displays HOMO energy, LUMO energy, and the HOMO-LUMO gap in both Hartrees (Eh) and electron Volts (eV).
    * **HOMO/LUMO Plotting:**
        * Generates a plot showing HOMO and LUMO energy levels (in eV) for all successfully parsed loaded files.
        * Allows customization of colors for the HOMO and LUMO series.
        * Uses `spiffyplot` for plot generation (output as SVG).

## Dependencies

* **Python 3.x**
* **Tkinter:** Usually included with standard Python installations.
* **`spiffyplot` (Recommended for Plotting):**
    * This script is configured to use `spiffyplot` for generating all plots.
    * `spiffyplot` is a command-line tool for creating scientific plots. It needs to be installed separately and be available in your system's PATH.
    * You can find `spiffyplot` here: [https://github.com/spiffyinput/spiffyplot](https://github.com/spiffyinput/spiffyplot)
* **Matplotlib (Fallback/Optional):**
    * The script checks for Matplotlib but primarily uses `spiffyplot`. If `spiffyplot` is not found, plotting features will be disabled. Matplotlib is not actively used for plotting in the current version if `spiffyplot` is the target.

## How to Run

1.  **Ensure Dependencies:**
    * Make sure Python 3 is installed.
    * Install `spiffyplot` and ensure it's in your system's PATH if you want plotting capabilities.
        ```bash
        # Example installation if using pip and spiffyplot is on PyPI (check spiffyplot docs)
        # pip install spiffyplot 
        # Or follow specific installation instructions from its repository.
        ```
2.  **Save the Script:** Save the Python code as a `.py` file (e.g., `orca_thermo_tools.py`).
3.  **Run from Terminal:**
    ```bash
    python orca_thermo_tools.py
    ```
    The GUI application window should appear.

## Using the Application

### Tab 1: Thermochemistry & Reaction Analysis

1.  **Load Files:** Click "Load ORCA File(s)" to select one or more ORCA output files (`.out`, `.log`) or single-frame `.xyz` files (for E\_el only).
    * Files will appear in the listbox. Files that could not be parsed or are missing essential data will be marked in red.
    * You can load additional files; they will be appended to the existing list.
2.  **Reorder Files:** Select a file in the list and use "Move Up ↑" or "Move Down ↓" to arrange them in the desired order for reaction profile plotting.
3.  **View Thermo Data:** Select a file from the list. Its detailed thermochemical data will appear in the main table on the right.
    * Use the "Energy Units" checkboxes to select which units (Eh, kcal/mol, kJ/mol, eV) are displayed.
    * Use the "Thermo Sections" checkboxes to toggle the visibility of different categories of thermochemical data.
4.  **Plot Reaction Profile:**
    * Select the "Plot Ref" (reference structure for relative energies) and "Plot Energy Type".
    * Choose a "Color" for the plot.
    * Click "Plot Profile". An SVG file will be generated by `spiffyplot` and the application will attempt to open it.
5.  **Curtin-Hammett Analysis:**
    * Select a "Reactant" from the dropdown (populated by loaded files).
    * Click "Add TS" to select one or more Transition State files. These will be added to the "Transition State(s)" listbox. You can remove or clear these.
    * Optionally, select a "Product" from the dropdown.
    * Choose a "Color" for the CH plot.
    * Click "Analyze Pathways" to see calculated $\Delta G^{\ddagger}$, $\Delta H^{\ddagger}$, etc., and product ratios in the table below.
    * Click "Plot CH Profile" to visualize the energy profile for the CH pathway using `spiffyplot`.
6.  **Step $\Delta$ Values:**
    * Click "Step $\Delta$ Values" to open a new window.
    * Select a "Reactant for Step" and "Product for Step" from the dropdowns.
    * Click "Add Step to Table". $\Delta G$, $\Delta H$, $\Delta S$, and $-T\Delta S$ for that step will be calculated and shown.
    * You can add multiple steps.
7.  **Consolidated Table:**
    * Click "Consolidated Table" to open a new window.
    * Select an "Energy Parameter" and "Display Unit".
    * The table will show the selected energy for all loaded files.

### Tab 2: TS Estimation (.allxyz)

1.  **Load .allxyz File(s):** Click "Load .allxyz / multi-frame .xyz File(s)" to select trajectory files.
2.  **Select File:** Choose a file from the "Loaded .allxyz Files" listbox.
3.  **View Frame Data:** The table on the right will show Frame #, Energy (Eh), and Relative Energy (kcal/mol) for each frame in the selected file. The highest energy frame is highlighted.
4.  **Save Highest E Geometry:** Click this button to save the coordinates of the highlighted highest energy frame as a new `.xyz` file.
5.  **Plot Profile:**
    * Choose a "Color" for the plot.
    * Click "Plot Profile for this .allxyz". An SVG file of Energy vs. Frame Number will be generated by `spiffyplot`.

### Tab 3: Orbital Analysis

1.  **Select File (in Tab 1):** The data displayed in this tab corresponds to the file currently selected in the "Loaded ORCA File(s)" list on Tab 1.
2.  **View Orbital Data:** The table shows HOMO energy, LUMO energy, and the HOMO-LUMO gap in Eh and eV for the selected file.
3.  **Plot HOMO/LUMO Levels:**
    * Choose colors for the HOMO and LUMO series.
    * Click "Plot HOMO/LUMO Levels". An SVG file showing HOMO and LUMO energies (eV) for all loaded files will be generated by `spiffyplot`.

## Notes on Spiffyplot

* This application relies on `spiffyplot` being installed and accessible in your system's PATH for all plotting functionalities.
* Plots are generated as SVG (Scalable Vector Graphics) files, typically in your system's temporary directory. The application will show a message with the path to the generated SVG and attempt to open it with your default SVG viewer (often a web browser).
* `spiffyplot` data files are created temporarily and then deleted.
* Line style and line width customizations for plots are simplified when using `spiffyplot` through this GUI, primarily focusing on color. `spiffyplot` itself has more advanced command-line options if needed.
* The Y-axis breaking feature, previously available with Matplotlib, is not directly supported by `spiffyplot` in the same way and has been removed from the UI for simplicity.

## Potential Issues

* **UnicodeDecodeError:** If you encounter this error when loading files, it means the file encoding is not standard UTF-8 or Latin-1. The script tries these two, but some ORCA outputs might have unusual characters.
* **Spiffyplot Not Found:** If `spiffyplot` is not installed or not in your PATH, plotting buttons will be disabled or will show an error.
* **Parsing Failures:** ORCA output formats can vary slightly. While the parsing logic tries to be robust, some files might not parse correctly, especially for orbital energies if they deviate significantly from common output patterns. Check the console for warnings.

This README provides a comprehensive guide to your application. Let me know if you'd like any adjustments or further details added!
