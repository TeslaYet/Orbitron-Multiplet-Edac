#!/usr/bin/env python3
"""
Multiplet2 GUI Application
This application provides a graphical interface for the Multiplet2 code,
allowing users to create input files, run the multiplet process, and
convert output files.
"""

import sys
import os
import subprocess
import re
import shutil
from PyQt6.QtWidgets import (QApplication, QMainWindow, QTabWidget, QWidget, QVBoxLayout,
                             QHBoxLayout, QLabel, QLineEdit, QPushButton, QTextEdit,
                             QFileDialog, QFormLayout, QGroupBox, QGridLayout, QMessageBox,
                             QSpinBox, QDoubleSpinBox, QScrollArea, QCheckBox, QTableWidget,
                             QTableWidgetItem, QComboBox)
from PyQt6.QtCore import Qt, QProcess, QThread, pyqtSignal
from PyQt6.QtGui import QPixmap

# Import the converter module
from convert_rpesalms import convert_rpesalms

# WorkerThread for running visualization commands
class WorkerThread(QThread):
    finished = pyqtSignal(bool, str)
    
    def __init__(self, command):
        super().__init__()
        self.command = command
        
    def run(self):
        try:
            output = subprocess.check_output(self.command, shell=True, stderr=subprocess.STDOUT)
            self.finished.emit(True, output.decode())
        except subprocess.CalledProcessError as e:
            self.finished.emit(False, e.output.decode())

class MultipletGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()
        
    def initUI(self):
        self.setWindowTitle('Orbitron-Multiplet-Edac')
        self.setGeometry(100, 100, 900, 700)
        
        # Main tab widget
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        
        # Create tabs
        self.input_tab = QWidget()
        self.run_tab = QWidget()
        self.convert_tab = QWidget()
        self.edac_tab = QWidget()
        self.viz_tab = QWidget()  # New visualization tab
        
        # Add tabs to widget
        self.tabs.addTab(self.input_tab, "Create Input")
        self.tabs.addTab(self.run_tab, "Run Multiplet")
        self.tabs.addTab(self.convert_tab, "Convert Output")
        self.tabs.addTab(self.edac_tab, "Run EDAC")
        self.tabs.addTab(self.viz_tab, "Visualization")  # Add new tab
        
        # Set up each tab
        self.setup_input_tab()
        self.setup_run_tab()
        self.setup_convert_tab()
        self.setup_edac_tab()
        self.setup_viz_tab()  # Setup the new visualization tab
    
    def setup_input_tab(self):
        # Main layout
        main_layout = QVBoxLayout()
        
        # Create scrollable text area to show the generated input
        self.input_preview = QTextEdit()
        self.input_preview.setReadOnly(True)
        
        # Create a form layout within a scroll area
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_content = QWidget()
        form_layout = QVBoxLayout(scroll_content)
        
        # Energy parameters group
        energy_box = QGroupBox("Energy Parameters")
        energy_layout = QGridLayout()
        
        # E(2p) and E(3d)
        self.e2p_input = QLineEdit("-639")
        self.e3d_input = QLineEdit("-1.e-6")
        energy_layout.addWidget(QLabel("E(2p):"), 0, 0)
        energy_layout.addWidget(self.e2p_input, 0, 1)
        energy_layout.addWidget(QLabel("E(3d):"), 0, 2)
        energy_layout.addWidget(self.e3d_input, 0, 3)
        
        # Crystal field matrix
        energy_layout.addWidget(QLabel("Crystal Field Matrix:"), 1, 0, 1, 4)
        self.cf_matrix = []
        for i in range(5):
            row = []
            for j in range(5):
                value = QLineEdit("0.")
                row.append(value)
                energy_layout.addWidget(value, i+2, j)
            self.cf_matrix.append(row)
        
        # Set crystal field values with exact spacing from reference
        self.cf_matrix[0][0].setText("-0.024")
        self.cf_matrix[0][1].setText("0.")
        self.cf_matrix[0][2].setText("0.")
        self.cf_matrix[0][3].setText("-0.056")
        self.cf_matrix[0][4].setText("0.")

        self.cf_matrix[1][0].setText("0.")
        self.cf_matrix[1][1].setText("0.064")
        self.cf_matrix[1][2].setText("0.")
        self.cf_matrix[1][3].setText("0.")
        self.cf_matrix[1][4].setText("0.056")

        self.cf_matrix[2][0].setText("0.")
        self.cf_matrix[2][1].setText("0.")
        self.cf_matrix[2][2].setText("-0.177")
        self.cf_matrix[2][3].setText("0.")
        self.cf_matrix[2][4].setText("0.")

        self.cf_matrix[3][0].setText("-0.056")
        self.cf_matrix[3][1].setText("0.")
        self.cf_matrix[3][2].setText("0.")
        self.cf_matrix[3][3].setText("0.064")
        self.cf_matrix[3][4].setText("0.")

        self.cf_matrix[4][0].setText("0.")
        self.cf_matrix[4][1].setText("0.056")
        self.cf_matrix[4][2].setText("0.")
        self.cf_matrix[4][3].setText("0.")
        self.cf_matrix[4][4].setText("-0.024")
        
        energy_box.setLayout(energy_layout)
        form_layout.addWidget(energy_box)
        
        # B-field parameters
        bfield_box = QGroupBox("Magnetic Field")
        bfield_layout = QHBoxLayout()
        self.bfield_strength = QLineEdit("1.e-3")
        self.bfield_theta = QLineEdit("90.")
        bfield_layout.addWidget(QLabel("B-field (eV):"))
        bfield_layout.addWidget(self.bfield_strength)
        bfield_layout.addWidget(QLabel("theta (deg):"))
        bfield_layout.addWidget(self.bfield_theta)
        bfield_box.setLayout(bfield_layout)
        form_layout.addWidget(bfield_box)
        
        # Photon energy parameters
        photon_box = QGroupBox("Photon Energy Settings")
        photon_layout = QHBoxLayout()
        self.omega_start = QLineEdit("651.8")
        self.omega_stop = QLineEdit("651.8")
        self.delta_omega = QLineEdit("2.")
        self.gamma = QLineEdit("0.4")
        self.gamma_flag = QLineEdit("0")
        photon_layout.addWidget(QLabel("ω start:"))
        photon_layout.addWidget(self.omega_start)
        photon_layout.addWidget(QLabel("ω stop:"))
        photon_layout.addWidget(self.omega_stop)
        photon_layout.addWidget(QLabel("Δω:"))
        photon_layout.addWidget(self.delta_omega)
        photon_layout.addWidget(QLabel("Γ:"))
        photon_layout.addWidget(self.gamma)
        photon_layout.addWidget(QLabel("flag:"))
        photon_layout.addWidget(self.gamma_flag)
        photon_box.setLayout(photon_layout)
        form_layout.addWidget(photon_box)
        
        # Electron configuration
        config_box = QGroupBox("Electron Configuration")
        config_layout = QVBoxLayout()
        
        config_row1 = QHBoxLayout()
        self.d_electrons = QLineEdit("5")
        config_row1.addWidget(QLabel("Number of d-electrons:"))
        config_row1.addWidget(self.d_electrons)
        config_layout.addLayout(config_row1)
        
        config_row2 = QHBoxLayout()
        self.l_values = QLineEdit("1 2 1 3 5")
        config_row2.addWidget(QLabel("l-values of active shells:"))
        config_row2.addWidget(self.l_values)
        config_layout.addLayout(config_row2)
        
        # SOC parameters
        config_row3 = QHBoxLayout()
        self.soc_params = QLineEdit("6.846  0.040 0 0 0")
        config_row3.addWidget(QLabel("SOC parameters:"))
        config_row3.addWidget(self.soc_params)
        config_layout.addLayout(config_row3)
        
        # Dipole matrix elements
        config_row4 = QHBoxLayout()
        self.dipole_elements = QLineEdit("2.064  0.02161  0.09695")
        config_row4.addWidget(QLabel("Dipole matrix elements:"))
        config_row4.addWidget(self.dipole_elements)
        config_layout.addLayout(config_row4)
        
        config_box.setLayout(config_layout)
        form_layout.addWidget(config_box)
        
        # Ground state parameters
        gs_box = QGroupBox("Ground State Parameters")
        gs_layout = QVBoxLayout()
        
        gs_row1 = QHBoxLayout()
        self.gs_config_count = QLineEdit("1")
        gs_row1.addWidget(QLabel("Number of configurations:"))
        gs_row1.addWidget(self.gs_config_count)
        gs_layout.addLayout(gs_row1)
        
        gs_row2 = QHBoxLayout()
        self.gs_occupation = QLineEdit("6 5 0 0 0")
        gs_row2.addWidget(QLabel("Occupation numbers:"))
        gs_row2.addWidget(self.gs_occupation)
        gs_layout.addLayout(gs_row2)
        
        # Slater-Condon parameters for ground state
        gs_row3 = QHBoxLayout()
        self.gs_slater_f2p3d = QLineEdit("0 0 5.0568")
        gs_row3.addWidget(QLabel("F_k(2p,3d):"))
        gs_row3.addWidget(self.gs_slater_f2p3d)
        gs_layout.addLayout(gs_row3)
        
        gs_row4 = QHBoxLayout()
        self.gs_slater_g2p3d = QLineEdit("0 3.6848 0 2.0936")
        gs_row4.addWidget(QLabel("G_k(2p,3d):"))
        gs_row4.addWidget(self.gs_slater_g2p3d)
        gs_layout.addLayout(gs_row4)
        
        gs_row5 = QHBoxLayout()
        self.gs_slater_f2p3d_2 = QLineEdit("0 0 5.0568")
        gs_row5.addWidget(QLabel("F_k(2p,3d) (2):"))
        gs_row5.addWidget(self.gs_slater_f2p3d_2)
        gs_layout.addLayout(gs_row5)
        
        gs_row6 = QHBoxLayout()
        self.gs_slater_f2p3d_3 = QLineEdit("0 0 5.0568")
        gs_row6.addWidget(QLabel("F_k(2p,3d) (3):"))
        gs_row6.addWidget(self.gs_slater_f2p3d_3)
        gs_layout.addLayout(gs_row6)
        
        gs_row7 = QHBoxLayout()
        self.gs_slater_g2p3d_2 = QLineEdit("0 3.6848 0 2.0936")
        gs_row7.addWidget(QLabel("G_k(2p,3d) (2):"))
        gs_row7.addWidget(self.gs_slater_g2p3d_2)
        gs_layout.addLayout(gs_row7)
        
        gs_row8 = QHBoxLayout()
        self.gs_slater_f3d3d = QLineEdit("0 0 9.4752 0 5.9256")
        gs_row8.addWidget(QLabel("F_k(3d,3d):"))
        gs_row8.addWidget(self.gs_slater_f3d3d)
        gs_layout.addLayout(gs_row8)
        
        gs_box.setLayout(gs_layout)
        form_layout.addWidget(gs_box)
        
        # Final state parameters
        fs_box = QGroupBox("Final State Parameters")
        fs_layout = QVBoxLayout()
        
        fs_row1 = QHBoxLayout()
        self.fs_config_count = QLineEdit("1")
        fs_row1.addWidget(QLabel("Number of configurations:"))
        fs_row1.addWidget(self.fs_config_count)
        fs_layout.addLayout(fs_row1)
        
        fs_row2 = QHBoxLayout()
        self.fs_occupation = QLineEdit("6 4 0 0 0")
        fs_row2.addWidget(QLabel("Occupation numbers:"))
        fs_row2.addWidget(self.fs_occupation)
        fs_layout.addLayout(fs_row2)
        
        # Slater-Condon parameters for final state
        fs_row3 = QHBoxLayout()
        self.fs_slater_f2p3d = QLineEdit("0 0 5.0568")
        fs_row3.addWidget(QLabel("F_k(2p,3d):"))
        fs_row3.addWidget(self.fs_slater_f2p3d)
        fs_layout.addLayout(fs_row3)
        
        fs_row4 = QHBoxLayout()
        self.fs_slater_g2p3d = QLineEdit("0 3.6848 0 2.0936")
        fs_row4.addWidget(QLabel("G_k(2p,3d):"))
        fs_row4.addWidget(self.fs_slater_g2p3d)
        fs_layout.addLayout(fs_row4)
        
        fs_row5 = QHBoxLayout()
        self.fs_slater_f2p3d_2 = QLineEdit("0 0 5.0568")
        fs_row5.addWidget(QLabel("F_k(2p,3d) (2):"))
        fs_row5.addWidget(self.fs_slater_f2p3d_2)
        fs_layout.addLayout(fs_row5)
        
        fs_row6 = QHBoxLayout()
        self.fs_slater_f2p3d_3 = QLineEdit("0 0 5.0568")
        fs_row6.addWidget(QLabel("F_k(2p,3d) (3):"))
        fs_row6.addWidget(self.fs_slater_f2p3d_3)
        fs_layout.addLayout(fs_row6)
        
        fs_row7 = QHBoxLayout()
        self.fs_slater_g2p3d_2 = QLineEdit("0 3.6848 0 2.0936")
        fs_row7.addWidget(QLabel("G_k(2p,3d) (2):"))
        fs_row7.addWidget(self.fs_slater_g2p3d_2)
        fs_layout.addLayout(fs_row7)
        
        fs_row8 = QHBoxLayout()
        self.fs_slater_f3d3d = QLineEdit("0 0 9.4752 0 5.9256")
        fs_row8.addWidget(QLabel("F_k(3d,3d):"))
        fs_row8.addWidget(self.fs_slater_f3d3d)
        fs_layout.addLayout(fs_row8)
        
        fs_box.setLayout(fs_layout)
        form_layout.addWidget(fs_box)
        
        # Intermediate state parameters
        is_box = QGroupBox("Intermediate State Parameters")
        is_layout = QVBoxLayout()
        
        is_row1 = QHBoxLayout()
        self.is_config_count = QLineEdit("1")
        is_row1.addWidget(QLabel("Number of configurations:"))
        is_row1.addWidget(self.is_config_count)
        is_layout.addLayout(is_row1)
        
        is_row2 = QHBoxLayout()
        self.is_occupation = QLineEdit("5 6  0 0 0")
        is_row2.addWidget(QLabel("Occupation numbers:"))
        is_row2.addWidget(self.is_occupation)
        is_layout.addLayout(is_row2)
        
        # Slater-Condon parameters for intermediate state
        is_row3 = QHBoxLayout()
        self.is_slater_f2p3d = QLineEdit("0 0 5.0568")
        is_row3.addWidget(QLabel("F_k(2p,3d):"))
        is_row3.addWidget(self.is_slater_f2p3d)
        is_layout.addLayout(is_row3)
        
        is_row4 = QHBoxLayout()
        self.is_slater_g2p3d = QLineEdit("0 3.6848 0 2.0936")
        is_row4.addWidget(QLabel("G_k(2p,3d):"))
        is_row4.addWidget(self.is_slater_g2p3d)
        is_layout.addLayout(is_row4)
        
        is_row5 = QHBoxLayout()
        self.is_slater_f2p3d_2 = QLineEdit("0 0 5.0568")
        is_row5.addWidget(QLabel("F_k(2p,3d) (2):"))
        is_row5.addWidget(self.is_slater_f2p3d_2)
        is_layout.addLayout(is_row5)
        
        is_row6 = QHBoxLayout()
        self.is_slater_f2p3d_3 = QLineEdit("0 0 5.0568")
        is_row6.addWidget(QLabel("F_k(2p,3d) (3):"))
        is_row6.addWidget(self.is_slater_f2p3d_3)
        is_layout.addLayout(is_row6)
        
        is_row7 = QHBoxLayout()
        self.is_slater_g2p3d_2 = QLineEdit("0 3.6848 0 2.0936")
        is_row7.addWidget(QLabel("G_k(2p,3d) (2):"))
        is_row7.addWidget(self.is_slater_g2p3d_2)
        is_layout.addLayout(is_row7)
        
        is_row8 = QHBoxLayout()
        self.is_slater_f3d3d = QLineEdit("0 0 8.9240 0 5.5528")
        is_row8.addWidget(QLabel("F_k(3d,3d):"))
        is_row8.addWidget(self.is_slater_f3d3d)
        is_layout.addLayout(is_row8)
        
        is_box.setLayout(is_layout)
        form_layout.addWidget(is_box)
        
        # Auger decay parameters
        auger_box = QGroupBox("Auger Decay Parameters")
        auger_layout = QVBoxLayout()
        
        auger_row1 = QHBoxLayout()
        self.auger_r_2pnp = QLineEdit("0. -.19047  0. -.15644")
        auger_row1.addWidget(QLabel("R_k(2p,np;3d,3d):"))
        auger_row1.addWidget(self.auger_r_2pnp)
        auger_layout.addLayout(auger_row1)
        
        auger_row2 = QHBoxLayout()
        self.auger_r_2pnf = QLineEdit("0. 0.7079   0. 0.44937")
        auger_row2.addWidget(QLabel("R_k(2p,nf;3d,3d):"))
        auger_row2.addWidget(self.auger_r_2pnf)
        auger_layout.addLayout(auger_row2)
        
        auger_row3 = QHBoxLayout()
        self.auger_r_2pnh = QLineEdit("0. 0.       0. 0.28778")
        auger_row3.addWidget(QLabel("R_k(2p,nh;3d,3d):"))
        auger_row3.addWidget(self.auger_r_2pnh)
        auger_layout.addLayout(auger_row3)
        
        auger_box.setLayout(auger_layout)
        form_layout.addWidget(auger_box)
        
        # Set the scroll area widget
        scroll_area.setWidget(scroll_content)
        main_layout.addWidget(scroll_area, 1)
        
        # Preview and save buttons
        button_layout = QHBoxLayout()
        self.preview_button = QPushButton("Preview Input")
        self.save_button = QPushButton("Save Input File")
        self.preview_button.clicked.connect(self.preview_input)
        self.save_button.clicked.connect(self.save_input_file)
        button_layout.addWidget(self.preview_button)
        button_layout.addWidget(self.save_button)
        main_layout.addLayout(button_layout)
        
        # Preview area
        preview_box = QGroupBox("Input File Preview")
        preview_layout = QVBoxLayout()
        preview_layout.addWidget(self.input_preview)
        preview_box.setLayout(preview_layout)
        main_layout.addWidget(preview_box)
        
        self.input_tab.setLayout(main_layout)
    
    def setup_run_tab(self):
        layout = QVBoxLayout()
        
        # Input file selection
        input_group = QGroupBox("Input File")
        input_layout = QHBoxLayout()
        self.input_file_path = QLineEdit()
        self.input_file_path.setReadOnly(True)
        browse_button = QPushButton("Browse")
        browse_button.clicked.connect(self.browse_input_file)
        input_layout.addWidget(QLabel("Input File:"))
        input_layout.addWidget(self.input_file_path)
        input_layout.addWidget(browse_button)
        input_group.setLayout(input_layout)
        layout.addWidget(input_group)
        
        # Output directory selection
        output_group = QGroupBox("Output Directory")
        output_layout = QHBoxLayout()
        self.output_dir_path = QLineEdit()
        self.output_dir_path.setReadOnly(True)
        output_browse_button = QPushButton("Browse")
        output_browse_button.clicked.connect(self.browse_output_dir)
        output_layout.addWidget(QLabel("Output Directory:"))
        output_layout.addWidget(self.output_dir_path)
        output_layout.addWidget(output_browse_button)
        output_group.setLayout(output_layout)
        layout.addWidget(output_group)
        
        # Run button
        self.run_button = QPushButton("Run Multiplet")
        self.run_button.clicked.connect(self.run_multiplet)
        layout.addWidget(self.run_button)
        
        # Console output
        console_group = QGroupBox("Console Output")
        console_layout = QVBoxLayout()
        self.console_output = QTextEdit()
        self.console_output.setReadOnly(True)
        console_layout.addWidget(self.console_output)
        console_group.setLayout(console_layout)
        layout.addWidget(console_group)
        
        self.run_tab.setLayout(layout)
        
        # Process for running the multiplet calculation
        self.process = QProcess()
        self.process.readyReadStandardOutput.connect(self.handle_stdout)
        self.process.readyReadStandardError.connect(self.handle_stderr)
        self.process.finished.connect(self.process_finished)
    
    def setup_convert_tab(self):
        layout = QVBoxLayout()
        
        # Input file selection
        input_group = QGroupBox("Input File (rpesalms.dat)")
        input_layout = QHBoxLayout()
        self.convert_input_path = QLineEdit()
        self.convert_input_path.setReadOnly(True)
        convert_browse_button = QPushButton("Browse")
        convert_browse_button.clicked.connect(self.browse_convert_input)
        input_layout.addWidget(QLabel("Input File:"))
        input_layout.addWidget(self.convert_input_path)
        input_layout.addWidget(convert_browse_button)
        input_group.setLayout(input_layout)
        layout.addWidget(input_group)
        
        # Output file selection
        output_group = QGroupBox("Output File (rpesalms.edac)")
        output_layout = QHBoxLayout()
        self.convert_output_path = QLineEdit()
        self.convert_output_path.setReadOnly(True)
        convert_output_button = QPushButton("Browse")
        convert_output_button.clicked.connect(self.browse_convert_output)
        output_layout.addWidget(QLabel("Output File:"))
        output_layout.addWidget(self.convert_output_path)
        output_layout.addWidget(convert_output_button)
        output_group.setLayout(output_layout)
        layout.addWidget(output_group)
        
        # Convert button
        self.convert_button = QPushButton("Convert File")
        self.convert_button.clicked.connect(self.convert_file)
        layout.addWidget(self.convert_button)
        
        # Status message
        self.convert_status = QLabel("Ready")
        layout.addWidget(self.convert_status)
        
        self.convert_tab.setLayout(layout)
    
    def setup_edac_tab(self):
        layout = QVBoxLayout()
        
        # Section 1: Select the existing rpesalms.edac file
        edac_files_group = QGroupBox("EDAC Files")
        edac_files_layout = QGridLayout()
        
        # File selection for rpesalms.edac
        self.edac_rpesalms_edac_path = QLineEdit()
        self.edac_rpesalms_edac_path.setReadOnly(True)
        rpesalms_edac_browse_button = QPushButton("Browse")
        rpesalms_edac_browse_button.clicked.connect(self.browse_edac_rpesalms_edac)
        edac_files_layout.addWidget(QLabel("rpesalms.edac File:"), 0, 0)
        edac_files_layout.addWidget(self.edac_rpesalms_edac_path, 0, 1)
        edac_files_layout.addWidget(rpesalms_edac_browse_button, 0, 2)
        
        # Cluster file selection
        self.edac_cluster_path = QLineEdit()
        self.edac_cluster_path.setReadOnly(True)
        cluster_browse_button = QPushButton("Browse")
        cluster_browse_button.clicked.connect(self.browse_edac_cluster)
        edac_files_layout.addWidget(QLabel("Cluster File (.clus):"), 1, 0)
        edac_files_layout.addWidget(self.edac_cluster_path, 1, 1)
        edac_files_layout.addWidget(cluster_browse_button, 1, 2)
        
        edac_files_group.setLayout(edac_files_layout)
        layout.addWidget(edac_files_group)
        
        # Section 2: Emitter Configuration
        emitter_group = QGroupBox("Emitter Configuration")
        emitter_layout = QVBoxLayout()
        
        # Number of emitters
        num_emitters_layout = QHBoxLayout()
        self.num_emitters_spin = QSpinBox()
        self.num_emitters_spin.setRange(1, 20)
        self.num_emitters_spin.setValue(4)  # Default value
        self.num_emitters_spin.valueChanged.connect(self.update_edac_emitter_table)
        
        num_emitters_layout.addWidget(QLabel("Number of emitters:"))
        num_emitters_layout.addWidget(self.num_emitters_spin)
        num_emitters_layout.addStretch()
        
        emitter_layout.addLayout(num_emitters_layout)
        
        # Emitter table
        self.emitter_table = QTableWidget(4, 2)  # Default 4 emitters
        self.emitter_table.setHorizontalHeaderLabels(["Emitter Number", "Emitter Type"])
        self.emitter_table.horizontalHeader().setStretchLastSection(True)
        
        # Set default values (from the original edac.in)
        default_values = ["1 25", "51 25", "296 25", "430 25"]
        for i, val in enumerate(default_values):
            parts = val.split()
            if len(parts) >= 2:
                self.emitter_table.setItem(i, 0, QTableWidgetItem(parts[0]))
                self.emitter_table.setItem(i, 1, QTableWidgetItem(parts[1]))
        
        emitter_layout.addWidget(self.emitter_table)
        emitter_group.setLayout(emitter_layout)
        
        layout.addWidget(emitter_group)
        
        # Section 3: EDAC Directory
        edac_dir_group = QGroupBox("EDAC Directory")
        edac_dir_layout = QHBoxLayout()
        self.edac_dir_path = QLineEdit()
        self.edac_dir_path.setReadOnly(True)
        
        # Default to "../Edac 2" relative to current file
        default_edac_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../Edac 2"))
        if os.path.exists(default_edac_path):
            self.edac_dir_path.setText(default_edac_path)
        
        edac_dir_browse_button = QPushButton("Browse")
        edac_dir_browse_button.clicked.connect(self.browse_edac_dir)
        edac_dir_layout.addWidget(QLabel("EDAC Directory:"))
        edac_dir_layout.addWidget(self.edac_dir_path)
        edac_dir_layout.addWidget(edac_dir_browse_button)
        edac_dir_group.setLayout(edac_dir_layout)
        layout.addWidget(edac_dir_group)
        
        # Run button
        self.run_edac_button = QPushButton("Run EDAC")
        self.run_edac_button.clicked.connect(self.run_edac)
        layout.addWidget(self.run_edac_button)
        
        # Console output
        console_group = QGroupBox("Console Output")
        console_layout = QVBoxLayout()
        self.edac_console_output = QTextEdit()
        self.edac_console_output.setReadOnly(True)
        console_layout.addWidget(self.edac_console_output)
        console_group.setLayout(console_layout)
        layout.addWidget(console_group)
        
        self.edac_tab.setLayout(layout)
        
        # Process for running EDAC commands
        self.edac_process = QProcess()
        self.edac_process.readyReadStandardOutput.connect(self.handle_edac_stdout)
        self.edac_process.readyReadStandardError.connect(self.handle_edac_stderr)
        self.edac_process.finished.connect(self.edac_process_finished)
    
    def generate_input_content(self):
        """Generate the content for the multiplet_input.txt file with exact formatting"""
        lines = []
        
        # Energy parameters - use actual input values while maintaining formatting
        lines.append(f"{self.e2p_input.text()} {self.e3d_input.text()}")
        
        # Crystal field matrix - use actual input values
        for i in range(5):
            row_values = []
            for j in range(5):
                # Add proper spacing to maintain format
                value = self.cf_matrix[i][j].text()
                # Pad with spaces to maintain alignment
                if len(value) < 7:
                    value = value.ljust(7)
                row_values.append(value)
            lines.append(" ".join(row_values))
        
        # B-field - use actual input values
        lines.append(f" {self.bfield_strength.text()} {self.bfield_theta.text()}")
        
        # Photon energy - use actual input values
        lines.append(f"{self.omega_start.text()}  {self.omega_stop.text()}   {self.delta_omega.text()} {self.gamma.text()} {self.gamma_flag.text()}")
        
        # Number of d-electrons
        lines.append(f"{self.d_electrons.text()}")
        
        # l-values - use actual input values
        lines.append(f"{self.l_values.text()}")
        
        # SOC parameters - use actual input values
        lines.append(f"{self.soc_params.text()}")
        
        # Dipole matrix elements - use actual input values
        lines.append(f"{self.dipole_elements.text()}")
        
        # Ground state configuration count
        lines.append(f"{self.gs_config_count.text()}")
        
        # Ground state occupation numbers
        lines.append(f"{self.gs_occupation.text()}")
        
        # Ground state Slater-Condon parameters - use actual input values
        lines.append(f"{self.gs_slater_f2p3d.text()}")
        lines.append(f"{self.gs_slater_g2p3d.text()}")
        lines.append(f"{self.gs_slater_f2p3d_2.text()}")
        lines.append(f"{self.gs_slater_f2p3d_3.text()}")
        lines.append(f"{self.gs_slater_g2p3d_2.text()}")
        lines.append(f"{self.gs_slater_f3d3d.text()}")
        
        # Final state configuration count
        lines.append(f"{self.fs_config_count.text()}")
        
        # Final state occupation numbers
        lines.append(f"{self.fs_occupation.text()}")
        
        # Final state Slater-Condon parameters - use actual input values
        lines.append(f"{self.fs_slater_f2p3d.text()}")
        lines.append(f"{self.fs_slater_g2p3d.text()}")
        lines.append(f"{self.fs_slater_f2p3d_2.text()}")
        lines.append(f"{self.fs_slater_f2p3d_3.text()}")
        lines.append(f"{self.fs_slater_g2p3d_2.text()}")
        lines.append(f"{self.fs_slater_f3d3d.text()}")
        
        # Intermediate state configuration count
        lines.append(f"{self.is_config_count.text()}")
        
        # Intermediate state occupation numbers
        lines.append(f"{self.is_occupation.text()}")
        
        # Intermediate state Slater-Condon parameters - use actual input values
        lines.append(f"{self.is_slater_f2p3d.text()}")
        lines.append(f"{self.is_slater_g2p3d.text()}")
        lines.append(f"{self.is_slater_f2p3d_2.text()}")
        lines.append(f"{self.is_slater_f2p3d_3.text()}")
        lines.append(f"{self.is_slater_g2p3d_2.text()}")
        lines.append(f"{self.is_slater_f3d3d.text()}")
        
        # Auger decay integrals - use actual input values
        lines.append(f"{self.auger_r_2pnp.text()}")
        lines.append(f"{self.auger_r_2pnf.text()}")
        lines.append(f"{self.auger_r_2pnh.text()}")
        
        return "\n".join(lines)
    
    def preview_input(self):
        """Preview the input file content"""
        content = self.generate_input_content()
        self.input_preview.setText(content)
    
    def save_input_file(self):
        """Save the input file with exact reference formatting"""
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Input File", 
                                                 "multiplet_input.txt", 
                                                 "Text Files (*.txt)")
        if file_path:
            content = self.generate_input_content()
            with open(file_path, 'w') as f:
                f.write(content)
            QMessageBox.information(self, "Success", 
                                   f"Input file saved to {file_path}\n\n"
                                   "Note: The file uses the exact reference values and formatting "
                                   "to ensure compatibility with the Multiplet code.")
    
    def browse_input_file(self):
        """Browse for input file"""
        file_path, _ = QFileDialog.getOpenFileName(self, "Select Input File", 
                                                  "", 
                                                  "Text Files (*.txt)")
        if file_path:
            self.input_file_path.setText(file_path)
    
    def browse_output_dir(self):
        """Browse for output directory"""
        dir_path = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir_path:
            self.output_dir_path.setText(dir_path)
    
    def run_multiplet(self):
        """Run the multiplet calculation"""
        input_file = self.input_file_path.text()
        output_dir = self.output_dir_path.text()
        
        if not input_file:
            QMessageBox.warning(self, "Error", "Please select an input file")
            return
            
        if not output_dir:
            QMessageBox.warning(self, "Error", "Please select an output directory")
            return
        
        # Make sure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Clear console output
        self.console_output.clear()
        
        # Disable run button during execution
        self.run_button.setEnabled(False)
        
        # Get path to multiplet executable (assumed to be in same directory as this script)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        multiplet_path = os.path.join(script_dir, "multiplet")
        
        # Change to output directory
        os.chdir(output_dir)
        
        # Start the process
        self.process.start(multiplet_path, [])
        
        # Send input file content to process stdin
        with open(input_file, 'r') as f:
            input_content = f.read()
        
        self.process.write(input_content.encode())
        self.process.closeWriteChannel()
    
    def handle_stdout(self):
        """Handle standard output from the process"""
        data = self.process.readAllStandardOutput().data().decode()
        self.console_output.append(data)
    
    def handle_stderr(self):
        """Handle standard error from the process"""
        data = self.process.readAllStandardError().data().decode()
        self.console_output.append(data)
    
    def process_finished(self, exit_code, exit_status):
        """Handle process completion"""
        self.run_button.setEnabled(True)
        
        if exit_code == 0:
            self.console_output.append("\nMultiplet calculation completed successfully!")
        else:
            self.console_output.append(f"\nMultiplet calculation failed with exit code {exit_code}")
    
    def browse_convert_input(self):
        """Browse for rpesalms.dat file to convert"""
        file_path, _ = QFileDialog.getOpenFileName(self, "Select rpesalms.dat File", 
                                                  "", 
                                                  "DAT Files (*.dat)")
        if file_path:
            self.convert_input_path.setText(file_path)
            
            # Suggest default output file path
            base_dir = os.path.dirname(file_path)
            base_name = os.path.basename(file_path)
            name_parts = os.path.splitext(base_name)
            default_output = os.path.join(base_dir, f"{name_parts[0]}.edac")
            self.convert_output_path.setText(default_output)
    
    def browse_convert_output(self):
        """Browse for rpesalms.edac output file location"""
        file_path, _ = QFileDialog.getSaveFileName(self, "Save Output File", 
                                                  "rpesalms.edac", 
                                                  "EDAC Files (*.edac)")
        if file_path:
            self.convert_output_path.setText(file_path)
    
    def convert_file(self):
        """Convert rpesalms.dat to rpesalms.edac"""
        input_file = self.convert_input_path.text()
        output_file = self.convert_output_path.text()
        
        if not input_file:
            QMessageBox.warning(self, "Error", "Please select an input file")
            return
            
        if not output_file:
            QMessageBox.warning(self, "Error", "Please specify an output file")
            return
        
        try:
            convert_rpesalms(input_file, output_file)
            self.convert_status.setText(f"Conversion successful! Output written to {output_file}")
            QMessageBox.information(self, "Success", f"File converted successfully!")
        except Exception as e:
            self.convert_status.setText(f"Conversion failed: {str(e)}")
            QMessageBox.critical(self, "Error", f"Conversion failed: {str(e)}")
    
    def browse_edac_rpesalms_edac(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, 
            "Select rpesalms.edac File", 
            os.path.dirname(os.path.abspath(__file__)),
            "EDAC Files (*.edac)"
        )
        if file_path:
            self.edac_rpesalms_edac_path.setText(file_path)
            self.check_edac_run_button()
    
    def browse_edac_cluster(self):
        file_path, _ = QFileDialog.getOpenFileName(
            self, 
            "Select Cluster File", 
            self.edac_dir_path.text() if self.edac_dir_path.text() else os.path.dirname(os.path.abspath(__file__)),
            "Cluster Files (*.clus)"
        )
        if file_path:
            self.edac_cluster_path.setText(file_path)
            self.check_edac_run_button()
    
    def browse_edac_dir(self):
        dir_path = QFileDialog.getExistingDirectory(
            self, 
            "Select EDAC Directory",
            self.edac_dir_path.text() if self.edac_dir_path.text() else os.path.dirname(os.path.abspath(__file__))
        )
        if dir_path:
            self.edac_dir_path.setText(dir_path)
            self.check_edac_run_button()
    
    def update_edac_emitter_table(self):
        num_emitters = self.num_emitters_spin.value()
        current_rows = self.emitter_table.rowCount()
        
        if num_emitters > current_rows:
            # Add new rows with default values
            for i in range(current_rows, num_emitters):
                self.emitter_table.insertRow(i)
                self.emitter_table.setItem(i, 0, QTableWidgetItem("1"))
                self.emitter_table.setItem(i, 1, QTableWidgetItem("25"))
        elif num_emitters < current_rows:
            # Remove excess rows
            for i in range(current_rows - 1, num_emitters - 1, -1):
                self.emitter_table.removeRow(i)
        
        self.check_edac_run_button()
    
    def check_edac_run_button(self):
        # Enable run button only if all required files are selected
        has_rpesalms = False
        has_cluster = False 
        has_dir = False
        
        if self.edac_rpesalms_edac_path.text() and os.path.exists(self.edac_rpesalms_edac_path.text()):
            has_rpesalms = True
        
        if self.edac_cluster_path.text() and os.path.exists(self.edac_cluster_path.text()):
            has_cluster = True
        
        if self.edac_dir_path.text() and os.path.exists(self.edac_dir_path.text()):
            has_dir = True
        
        self.run_edac_button.setEnabled(has_rpesalms and has_cluster and has_dir)
    
    def run_edac(self):
        # Validate inputs
        if not self.edac_rpesalms_edac_path.text() or not os.path.exists(self.edac_rpesalms_edac_path.text()):
            QMessageBox.warning(self, "Missing Input", "Please select a rpesalms.edac file")
            return
        
        if not self.edac_cluster_path.text() or not os.path.exists(self.edac_cluster_path.text()):
            QMessageBox.warning(self, "Missing Input", "Please select a cluster file")
            return
        
        if not self.edac_dir_path.text() or not os.path.exists(self.edac_dir_path.text()):
            QMessageBox.warning(self, "Invalid Path", "EDAC directory does not exist")
            return
        
        # Clear console and disable button during processing
        self.edac_console_output.clear()
        self.run_edac_button.setEnabled(False)
        self.edac_console_output.append("Starting EDAC workflow...\n")
        
        # Store the original directory to return to later
        original_dir = os.getcwd()
        
        try:
            # Change to EDAC directory for the entire process
            os.chdir(self.edac_dir_path.text())
            self.edac_console_output.append(f"Changed to directory: {os.getcwd()}")
            
            # Copy the selected rpesalms.edac file to the EDAC directory
            rpesalms_edac = self.edac_rpesalms_edac_path.text()
            edac_target = os.path.join(os.getcwd(), "rpesalms.edac")
            
            self.edac_console_output.append(f"Copying {rpesalms_edac} to {edac_target}")
            try:
                shutil.copy2(rpesalms_edac, edac_target)
                self.edac_console_output.append("File copied successfully")
            except Exception as e:
                self.edac_console_output.append(f"Error copying rpesalms.edac file: {e}")
                os.chdir(original_dir)  # Return to original directory
                self.run_edac_button.setEnabled(True)
                return
            
            # Copy the selected cluster file to the EDAC directory
            cluster_file = self.edac_cluster_path.text()
            cluster_filename = os.path.basename(cluster_file)
            cluster_target = os.path.join(os.getcwd(), cluster_filename)
            
            self.edac_console_output.append(f"Copying {cluster_file} to {cluster_target}")
            try:
                shutil.copy2(cluster_file, cluster_target)
                self.edac_console_output.append("File copied successfully")
            except Exception as e:
                self.edac_console_output.append(f"Error copying cluster file: {e}")
                os.chdir(original_dir)  # Return to original directory
                self.run_edac_button.setEnabled(True)
                return
            
            # Create rpededac input file - like in the EDAC GUI
            input_file = "rpededac_input.txt"
            with open(input_file, "w") as f:
                # Write cluster filename
                f.write(f"{cluster_filename}\n")
                
                # Write number of emitters
                num_emitters = self.num_emitters_spin.value()
                f.write(f"{num_emitters}\n")
                
                # Write emitter properties
                for i in range(num_emitters):
                    emitter_num = self.emitter_table.item(i, 0).text()
                    emitter_type = self.emitter_table.item(i, 1).text()
                    f.write(f"{emitter_num} {emitter_type} ")
                f.write("\n")
            
            self.edac_console_output.append(f"Created rpededac input file: {input_file}")
            
            # Run rpededac like in the EDAC GUI
            self.edac_console_output.append("\nRunning rpededac (which will run EDAC multiple times)...")
            
            # Platform-specific executable name
            if sys.platform == "darwin":  # macOS
                rpededac_cmd = "./rpededac"
            else:  # Linux or Windows
                # Check if rpededac.exe exists
                if os.path.exists("rpededac.exe"):
                    rpededac_cmd = "./rpededac.exe"
                else:
                    rpededac_cmd = "./rpededac"
            
            # Use cat to pipe the input file to rpededac - exactly like the EDAC GUI
            cmd = f"cat {input_file} | {rpededac_cmd}"
            self.edac_console_output.append(f"Executing command: {cmd}")
            
            self.edac_worker = WorkerThread(cmd)
            self.edac_worker.finished.connect(self.edac_worker_finished)
            self.edac_worker.start()
            
        except Exception as e:
            self.edac_console_output.append(f"Error in EDAC workflow: {e}")
            # Return to original directory
            os.chdir(original_dir)
            self.run_edac_button.setEnabled(True)
    
    def edac_worker_finished(self, success, output):
        """Handle completion of EDAC calculation from WorkerThread"""
        # This is called when the WorkerThread completes
        if success:
            self.edac_console_output.append("EDAC calculation completed successfully")
            self.edac_console_output.append("Output: " + output[:200] + "..." if len(output) > 200 else output)
            
            # Check for generated MS files
            ms_files = [f for f in os.listdir('.') if f.endswith('.ms')]
            if ms_files:
                self.edac_console_output.append(f"Generated {len(ms_files)} .ms files: {', '.join(ms_files[:5])}{'...' if len(ms_files) > 5 else ''}")
                
                # Update visualization tab settings
                # Keep current directory for visualization
                self.viz_dir_path.setText(os.getcwd())
                
                # Automatically switch to visualization tab
                self.tabs.setCurrentIndex(4)  # Switch to visualization tab (index 4)
                
                # Refresh MS files in viz tab
                self.refresh_ms_files()
                
                # Show success message
                QMessageBox.information(self, "Success", f"EDAC calculation completed successfully! Generated {len(ms_files)} .ms files.\nSwitched to Visualization tab.")
            else:
                self.edac_console_output.append("Warning: No .ms files were generated")
                QMessageBox.warning(self, "Warning", "EDAC calculation completed but no .ms files were generated.")
        else:
            self.edac_console_output.append(f"EDAC calculation failed with error: {output}")
            QMessageBox.critical(self, "Error", f"EDAC calculation failed:\n{output[:500]}...")
        
        # Re-enable the run button
        self.run_edac_button.setEnabled(True)
    
    def create_edac_input(self):
        """Create the edac.in file with parameters from the GUI"""
        cluster_filename = os.path.basename(self.edac_cluster_path.text())
        edac_in_path = os.path.join(self.edac_dir_path.text(), "edac.in")
        
        try:
            with open(edac_in_path, 'w') as f:
                # Cluster information
                f.write(f"cluster input {cluster_filename}\n")
                f.write("cluster surface on\n")
                
                # Default emission angle settings
                f.write("emission angle phi 0 360 181\n")
                f.write("emission angle theta 0 80 41\n")
                f.write("emission angle window 3\n")
                
                # Default energy settings
                f.write("emission energy E(eV) 622.4031 622.4031 1\n")
                
                # Emitters from the table
                emitter_str = "emitters " + str(self.num_emitters_spin.value())
                for i in range(self.num_emitters_spin.value()):
                    emitter_num = self.emitter_table.item(i, 0).text()
                    emitter_type = self.emitter_table.item(i, 1).text()
                    emitter_str += f" {emitter_num} {emitter_type}"
                f.write(f"{emitter_str}\n")
                
                # Default settings
                f.write("muffin-tin\n")
                f.write("orders 1 5\n")
                f.write("V0 E(eV) 10.\n")
                f.write("imfp imfp_input\n")
                
                # WF manual section with actual coefficients
                f.write("wf manual on\n")
                f.write("wf manual 2 5  0 0 7.20988e-06 2.57039e-06 8.68755e-06 -7.07164e-06 3.55705e-07 -8.98617e-06 0 0 0 0 0 0 0 0 0 0 -0.000129915 -0.000234869 0.000316588 0.000571304 -0.000489105 -0.000721222 0.000609963 0.000761132 -0.000561167 -0.000693942 0.000352114 0.000484381 -0.000109304 -0.000168153 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -7.79543e-07 1.51005e-07 2.12498e-08 1.47513e-06 2.90524e-06 -2.75591e-06 4.06886e-07 -1.52056e-06 -7.99225e-06 6.11839e-06 5.31244e-06 -5.0541e-06 -1.63129e-08 -9.44899e-08 -5.05541e-06 2.34583e-06 1.24108e-06 -1.44111e-06 -8.48157e-07 2.37151e-06 3.32398e-06 -4.76868e-06 0 0 4.25226e-06 3.27374e-06 -5.52058e-06 8.86954e-06 -1.64526e-05 7.5907e-07 0 0 0 0 0 0 0 0 0 0 0.000136687 0.000223885 -0.000370661 -0.000501411 0.00055359 0.000669528 -0.000604168 -0.000775558 0.000504584 0.00073476 -0.000335665 -0.000481788 0.000146973 0.000115325 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 3.81725e-06 -4.53819e-06 -2.3211e-06 5.11055e-06 -4.29208e-06 2.93305e-06 8.41215e-06 -8.51404e-06 -4.46301e-06 5.15462e-06 -3.45587e-06 9.7106e-07 8.70031e-06 -2.51351e-06 -1.75053e-06 1.28941e-06 7.25956e-07 -1.19412e-06 9.80157e-07 1.35474e-06 -7.62385e-07 1.30065e-07\n")
                
                # Default output settings
                f.write("report emitters\n")
                f.write("scan pd o_0_0_209.ms\n")
                f.write("report size\n")
                f.write("end\n")
            
            self.edac_console_output.append(f"Created edac.in file at {edac_in_path}")
            return True
        except Exception as e:
            self.edac_console_output.append(f"Error creating edac.in file: {e}")
            return False
    
    def handle_edac_stdout(self):
        data = self.edac_process.readAllStandardOutput().data().decode()
        self.edac_console_output.append(data)
    
    def handle_edac_stderr(self):
        data = self.edac_process.readAllStandardError().data().decode()
        self.edac_console_output.append(data)
    
    def edac_process_finished(self, exit_code, exit_status):
        if exit_code == 0:
            self.edac_console_output.append("EDAC calculation completed successfully")
            QMessageBox.information(self, "Success", "EDAC calculation completed successfully!")
        else:
            self.edac_console_output.append(f"EDAC calculation failed with exit code {exit_code}")
            QMessageBox.critical(self, "Error", f"EDAC calculation failed with exit code {exit_code}")
        
        self.run_edac_button.setEnabled(True)

    def setup_viz_tab(self):
        """Setup the visualization tab for viewing .ms files with intens_stereo tools"""
        layout = QVBoxLayout()
        
        # MS file selection
        ms_group = QGroupBox("MS File Selection")
        ms_layout = QGridLayout()
        
        self.ms_file_label = QLabel("MS file:")
        self.ms_file_combo = QComboBox()
        self.ms_file_combo.currentIndexChanged.connect(self.update_preview)
        self.refresh_ms_btn = QPushButton("Refresh")
        self.refresh_ms_btn.clicked.connect(self.refresh_ms_files)
        
        ms_layout.addWidget(self.ms_file_label, 0, 0)
        ms_layout.addWidget(self.ms_file_combo, 0, 1)
        ms_layout.addWidget(self.refresh_ms_btn, 0, 2)
        
        ms_group.setLayout(ms_layout)
        layout.addWidget(ms_group)
        
        # Directory selection
        dir_group = QGroupBox("Working Directory")
        dir_layout = QHBoxLayout()
        self.viz_dir_path = QLineEdit()
        self.viz_dir_path.setReadOnly(True)
        
        # Default to Edac 2 directory
        default_edac_path = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../Edac 2"))
        if os.path.exists(default_edac_path):
            self.viz_dir_path.setText(default_edac_path)
        
        dir_browse_button = QPushButton("Browse")
        dir_browse_button.clicked.connect(self.browse_viz_dir)
        dir_layout.addWidget(QLabel("Directory:"))
        dir_layout.addWidget(self.viz_dir_path)
        dir_layout.addWidget(dir_browse_button)
        dir_group.setLayout(dir_layout)
        layout.addWidget(dir_group)
        
        # Visualization options
        viz_group = QGroupBox("Visualization Options")
        viz_layout = QVBoxLayout()
        
        # Type selection
        type_layout = QHBoxLayout()
        self.viz_type_label = QLabel("Visualization type:")
        self.viz_type_combo = QComboBox()
        self.viz_type_combo.addItems(["Hot colormap (intens_stereo_hot)", 
                                      "Red-Blue colormap (intens_stereo_rb)"])
        
        type_layout.addWidget(self.viz_type_label)
        type_layout.addWidget(self.viz_type_combo)
        type_layout.addStretch()
        
        viz_layout.addLayout(type_layout)
        
        # Generate button
        self.generate_btn = QPushButton("Generate Visualization")
        self.generate_btn.clicked.connect(self.generate_visualization)
        self.generate_btn.setEnabled(False)
        
        viz_layout.addWidget(self.generate_btn)
        
        viz_group.setLayout(viz_layout)
        layout.addWidget(viz_group)
        
        # Status label
        self.viz_status_label = QLabel("Status: Ready")
        layout.addWidget(self.viz_status_label)
        
        # Image display
        image_group = QGroupBox("Image Preview")
        image_layout = QVBoxLayout()
        
        self.image_label = QLabel("No image generated yet")
        self.image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.scroll_area = QScrollArea()
        self.scroll_area.setWidget(self.image_label)
        self.scroll_area.setWidgetResizable(True)
        
        image_layout.addWidget(self.scroll_area)
        
        image_group.setLayout(image_layout)
        layout.addWidget(image_group)
        
        self.viz_tab.setLayout(layout)
    
    def browse_viz_dir(self):
        """Browse for the visualization working directory"""
        dir_path = QFileDialog.getExistingDirectory(
            self, 
            "Select Working Directory",
            self.viz_dir_path.text() if self.viz_dir_path.text() else os.path.dirname(os.path.abspath(__file__))
        )
        if dir_path:
            self.viz_dir_path.setText(dir_path)
            self.refresh_ms_files()
    
    def refresh_ms_files(self):
        """Refresh the list of .ms files in the current directory"""
        self.ms_file_combo.clear()
        
        if not self.viz_dir_path.text() or not os.path.exists(self.viz_dir_path.text()):
            return
            
        # Find all .ms files in the selected directory
        ms_files = [f for f in os.listdir(self.viz_dir_path.text()) if f.endswith('.ms')]
        
        if ms_files:
            self.ms_file_combo.addItems(ms_files)
            self.generate_btn.setEnabled(True)
        else:
            self.generate_btn.setEnabled(False)
    
    def update_preview(self):
        """Update the image preview when a new MS file is selected"""
        if not self.viz_dir_path.text() or not os.path.exists(self.viz_dir_path.text()):
            return
            
        # Check if there's already a PNG file for this MS file
        selected_file = self.ms_file_combo.currentText()
        if not selected_file:
            return
            
        # Check for existing PNG files with pattern matching the base name
        base_name = os.path.splitext(selected_file)[0]
        possible_files = [f for f in os.listdir(self.viz_dir_path.text()) 
                         if f.startswith(base_name) and f.endswith('.png')]
        
        if possible_files:
            # Use the most recent matching PNG file
            possible_files.sort(key=lambda f: os.path.getmtime(os.path.join(self.viz_dir_path.text(), f)), 
                             reverse=True)
            self.display_image(os.path.join(self.viz_dir_path.text(), possible_files[0]))
    
    def generate_visualization(self):
        """Generate visualization for the selected MS file"""
        selected_file = self.ms_file_combo.currentText()
        if not selected_file:
            return
        
        # Debug directory info
        self.viz_status_label.setText(f"Status: Working with directory: '{self.viz_dir_path.text()}'")
        
        # Verify working directory
        if not os.path.exists(self.viz_dir_path.text()):
            QMessageBox.critical(self, "Error", f"Working directory does not exist: '{self.viz_dir_path.text()}'")
            return
        
        # List files in directory to debug
        try:
            dir_files = os.listdir(self.viz_dir_path.text())
            executable_files = [f for f in dir_files if f.startswith("intens_stereo")]
            self.viz_status_label.setText(f"Status: Found visualization files: {executable_files}")
        except Exception as e:
            self.viz_status_label.setText(f"Status: Error listing directory: {e}")
            QMessageBox.critical(self, "Error", f"Error listing directory: {e}")
            return
        
        # Determine which visualization program to use
        if self.viz_type_combo.currentIndex() == 0:
            # Hot colormap
            viz_prog = "intens_stereo_hot"
        else:
            # Red-blue colormap
            viz_prog = "intens_stereo_rb"
        
        # Find the executable with explicit paths
        found_executable = False
        executable_path = ""
        
        # Check for the executable with exact name
        if viz_prog in executable_files:
            executable_path = os.path.join(self.viz_dir_path.text(), viz_prog)
            found_executable = True
        # Check for executable with .exe extension
        elif f"{viz_prog}.exe" in executable_files:
            executable_path = os.path.join(self.viz_dir_path.text(), f"{viz_prog}.exe")
            found_executable = True
            
        if not found_executable:
            error_msg = f"Could not find '{viz_prog}' or '{viz_prog}.exe' in '{self.viz_dir_path.text()}'. Available files: {executable_files}"
            self.viz_status_label.setText(f"Status: Error - {error_msg}")
            QMessageBox.critical(self, "Error", error_msg)
            return
        
        # Make sure the selected MS file exists in the working directory
        ms_file_path = os.path.join(self.viz_dir_path.text(), selected_file)
        if not os.path.exists(ms_file_path):
            error_msg = f"MS file not found: '{ms_file_path}'"
            self.viz_status_label.setText(f"Status: Error - {error_msg}")
            QMessageBox.critical(self, "Error", error_msg)
            return
        
        # Save current directory and change to working directory
        current_dir = os.getcwd()
        try:
            os.chdir(self.viz_dir_path.text())
            self.viz_status_label.setText(f"Status: Changed to directory: '{os.getcwd()}'")
            
            # Build command - use relative paths since we're in the correct directory
            if sys.platform == "darwin":  # macOS
                cmd = f"./{os.path.basename(executable_path)} {selected_file}"
            else:  # Linux
                cmd = f"./{os.path.basename(executable_path)} {selected_file}"
            
            self.viz_status_label.setText(f"Status: Running command: '{cmd}' in '{os.getcwd()}'")
            self.generate_btn.setEnabled(False)
            
            # Run the command
            self.worker = WorkerThread(cmd)
            self.worker.finished.connect(self.visualization_finished)
            self.worker.start()
            
        except Exception as e:
            error_msg = f"Error running visualization: {e}"
            self.viz_status_label.setText(f"Status: Error - {error_msg}")
            QMessageBox.critical(self, "Error", error_msg)
            # Return to original directory
            os.chdir(current_dir)
            self.generate_btn.setEnabled(True)
    
    def visualization_finished(self, success, output):
        """Handle completion of visualization generation"""
        # Return to original directory if needed
        if os.getcwd() != os.path.dirname(os.path.abspath(__file__)):
            try:
                os.chdir(os.path.dirname(os.path.abspath(__file__)))
            except Exception:
                pass

        if success:
            self.viz_status_label.setText(f"Status: Visualization generated successfully. Output: {output[:100]}...")
            selected_file = self.ms_file_combo.currentText()
            base_name = os.path.splitext(selected_file)[0]
            
            # Look for the newly created PNG file
            possible_files = [f for f in os.listdir(self.viz_dir_path.text()) 
                             if f.startswith(base_name) and f.endswith('.png')]
            
            if possible_files:
                # Use the most recent matching PNG file
                possible_files.sort(key=lambda f: os.path.getmtime(os.path.join(self.viz_dir_path.text(), f)), 
                                 reverse=True)
                self.display_image(os.path.join(self.viz_dir_path.text(), possible_files[0]))
            else:
                self.viz_status_label.setText(f"Status: Visualization completed but no image found. Files in dir: {os.listdir(self.viz_dir_path.text())[:10]}")
                QMessageBox.warning(self, "Warning", "Visualization generated but no output image found.")
        else:
            self.viz_status_label.setText(f"Status: Error generating visualization: {output}")
            QMessageBox.critical(self, "Error", f"Error generating visualization:\n{output}")
        
        self.generate_btn.setEnabled(True)
    
    def display_image(self, image_path):
        """Display an image in the visualization tab"""
        pixmap = QPixmap(image_path)
        if not pixmap.isNull():
            self.image_label.setPixmap(pixmap)
            self.image_label.resize(pixmap.size())
            self.image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        else:
            self.image_label.setText(f"Cannot load image: {image_path}")
            self.viz_status_label.setText(f"Status: Failed to load image: {image_path}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    ex = MultipletGUI()
    ex.show()
    sys.exit(app.exec()) 