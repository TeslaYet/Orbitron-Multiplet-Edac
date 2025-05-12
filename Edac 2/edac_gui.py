import sys
import os
import subprocess
from PyQt6.QtWidgets import (QApplication, QMainWindow, QTabWidget, QWidget, QLabel, 
                            QLineEdit, QPushButton, QFileDialog, QVBoxLayout, QHBoxLayout, 
                            QGridLayout, QSpinBox, QTableWidget, QTableWidgetItem, 
                            QScrollArea, QComboBox, QGroupBox, QMessageBox)
from PyQt6.QtGui import QPixmap, QImage
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QLibraryInfo

# Fix macOS plugin path issues
if sys.platform == "darwin":
    # Set the plugin path for macOS
    os.environ['QT_QPA_PLATFORM_PLUGIN_PATH'] = QLibraryInfo.path(QLibraryInfo.LibraryPath.PluginsPath)
    # Debug output for troubleshooting
    print(f"Plugin path set to: {os.environ['QT_QPA_PLATFORM_PLUGIN_PATH']}")

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

class EdacGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("EDAC Interface")
        self.setGeometry(100, 100, 800, 600)
        
        # Initialize variables
        self.rpesalms_path = ""
        self.cluster_path = ""
        self.output_dir = os.getcwd()
        self.output_files = []
        
        # Create tabs
        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)
        
        # Create the main input tab
        self.input_tab = QWidget()
        self.tabs.addTab(self.input_tab, "RPEDEDAC Input")
        
        # Create the visualization tab
        self.viz_tab = QWidget()
        self.tabs.addTab(self.viz_tab, "Visualization")
        
        self.setup_input_tab()
        self.setup_viz_tab()
        
    def setup_input_tab(self):
        layout = QVBoxLayout()
        
        # File selection section
        file_group = QGroupBox("Input Files")
        file_layout = QGridLayout()
        
        # rpesalms.edac selection
        self.rpesalms_label = QLabel("rpesalms.edac file:")
        self.rpesalms_path_edit = QLineEdit()
        self.rpesalms_path_edit.setReadOnly(True)
        self.rpesalms_browse_btn = QPushButton("Browse")
        self.rpesalms_browse_btn.clicked.connect(self.browse_rpesalms)
        
        file_layout.addWidget(self.rpesalms_label, 0, 0)
        file_layout.addWidget(self.rpesalms_path_edit, 0, 1)
        file_layout.addWidget(self.rpesalms_browse_btn, 0, 2)
        
        # Cluster file selection
        self.cluster_label = QLabel("Cluster file (.clus):")
        self.cluster_path_edit = QLineEdit()
        self.cluster_path_edit.setReadOnly(True)
        self.cluster_browse_btn = QPushButton("Browse")
        self.cluster_browse_btn.clicked.connect(self.browse_cluster)
        
        file_layout.addWidget(self.cluster_label, 1, 0)
        file_layout.addWidget(self.cluster_path_edit, 1, 1)
        file_layout.addWidget(self.cluster_browse_btn, 1, 2)
        
        file_group.setLayout(file_layout)
        layout.addWidget(file_group)
        
        # Emitter configuration
        emitter_group = QGroupBox("Emitter Configuration")
        emitter_layout = QVBoxLayout()
        
        # Number of emitters
        num_emitters_layout = QHBoxLayout()
        self.num_emitters_label = QLabel("Number of emitters:")
        self.num_emitters_spin = QSpinBox()
        self.num_emitters_spin.setRange(1, 20)
        self.num_emitters_spin.valueChanged.connect(self.update_emitter_table)
        
        num_emitters_layout.addWidget(self.num_emitters_label)
        num_emitters_layout.addWidget(self.num_emitters_spin)
        num_emitters_layout.addStretch()
        
        emitter_layout.addLayout(num_emitters_layout)
        
        # Emitter table
        self.emitter_table = QTableWidget(0, 2)
        self.emitter_table.setHorizontalHeaderLabels(["Emitter Number", "Emitter Type"])
        self.emitter_table.horizontalHeader().setStretchLastSection(True)
        
        emitter_layout.addWidget(self.emitter_table)
        emitter_group.setLayout(emitter_layout)
        
        layout.addWidget(emitter_group)
        
        # Run button
        self.run_btn = QPushButton("Run RPEDEDAC")
        self.run_btn.clicked.connect(self.run_rpededac)
        self.run_btn.setEnabled(False)
        
        layout.addWidget(self.run_btn)
        
        # Status label
        self.status_label = QLabel("Status: Ready")
        layout.addWidget(self.status_label)
        
        self.input_tab.setLayout(layout)
        
    def setup_viz_tab(self):
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
        
    def browse_rpesalms(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select rpesalms.edac file", "", "EDAC Files (*.edac)")
        if file_path:
            self.rpesalms_path = file_path
            self.rpesalms_path_edit.setText(file_path)
            self.check_run_button()
    
    def browse_cluster(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select cluster file", "", "Cluster Files (*.clus)")
        if file_path:
            self.cluster_path = file_path
            self.cluster_path_edit.setText(file_path)
            self.check_run_button()
    
    def update_emitter_table(self):
        num_emitters = self.num_emitters_spin.value()
        self.emitter_table.setRowCount(num_emitters)
        
        # Preserve existing values
        for i in range(num_emitters):
            # Create new items if they don't exist
            if self.emitter_table.item(i, 0) is None:
                self.emitter_table.setItem(i, 0, QTableWidgetItem("1"))
            if self.emitter_table.item(i, 1) is None:
                self.emitter_table.setItem(i, 1, QTableWidgetItem("25"))
        
        self.check_run_button()
    
    def check_run_button(self):
        if (self.rpesalms_path and self.cluster_path and self.num_emitters_spin.value() > 0):
            self.run_btn.setEnabled(True)
        else:
            self.run_btn.setEnabled(False)
    
    def run_rpededac(self):
        # Create input file for rpededac
        input_file = "rpededac_input.txt"
        with open(input_file, "w") as f:
            # Write cluster filename
            cluster_filename = os.path.basename(self.cluster_path)
            f.write(f"{cluster_filename}\n")
            
            # Write number of emitters
            num_emitters = self.num_emitters_spin.value()
            f.write(f"{num_emitters}\n")
            
            # Write emitter properties
            emitter_data = []
            for i in range(num_emitters):
                emitter_num = self.emitter_table.item(i, 0).text()
                emitter_type = self.emitter_table.item(i, 1).text()
                f.write(f"{emitter_num} {emitter_type} ")
            f.write("\n")
        
        # Copy necessary files to working directory if needed
        if os.path.dirname(self.rpesalms_path) != os.getcwd():
            cmd_copy_rpesalms = f"cp \"{self.rpesalms_path}\" rpesalms.edac"
            subprocess.run(cmd_copy_rpesalms, shell=True)
        
        if os.path.dirname(self.cluster_path) != os.getcwd():
            cmd_copy_cluster = f"cp \"{self.cluster_path}\" {cluster_filename}"
            subprocess.run(cmd_copy_cluster, shell=True)
        
        # Run rpededac with input
        self.status_label.setText("Status: Running RPEDEDAC...")
        self.run_btn.setEnabled(False)
        
        cmd = f"cat {input_file} | ./rpededac"
        self.worker = WorkerThread(cmd)
        self.worker.finished.connect(self.rpededac_finished)
        self.worker.start()
    
    def rpededac_finished(self, success, output):
        if success:
            self.status_label.setText("Status: RPEDEDAC finished successfully")
            self.refresh_ms_files()
            QMessageBox.information(self, "Success", "RPEDEDAC finished successfully!")
        else:
            self.status_label.setText("Status: Error running RPEDEDAC")
            QMessageBox.critical(self, "Error", f"Error running RPEDEDAC:\n{output}")
        
        self.run_btn.setEnabled(True)
    
    def refresh_ms_files(self):
        self.ms_file_combo.clear()
        
        # Find all .ms files in current directory
        ms_files = [f for f in os.listdir('.') if f.endswith('.ms')]
        self.output_files = ms_files
        
        if ms_files:
            self.ms_file_combo.addItems(ms_files)
            self.generate_btn.setEnabled(True)
        else:
            self.generate_btn.setEnabled(False)
    
    def update_preview(self):
        # If there's an existing visualization for this file, show it
        selected_file = self.ms_file_combo.currentText()
        if not selected_file:
            return
            
        # Check if there's already a PNG file for this MS file
        base_name = os.path.splitext(selected_file)[0]
        possible_files = [f for f in os.listdir('.') if f.startswith(base_name) and f.endswith('.png')]
        
        if possible_files:
            self.display_image(possible_files[0])
    
    def generate_visualization(self):
        selected_file = self.ms_file_combo.currentText()
        if not selected_file:
            return
        
        # Determine which visualization program to use
        if self.viz_type_combo.currentIndex() == 0:
            cmd = f"./intens_stereo_hot {selected_file}"
        else:
            cmd = f"./intens_stereo_rb {selected_file}"
        
        self.status_label.setText(f"Status: Generating visualization for {selected_file}...")
        self.generate_btn.setEnabled(False)
        
        self.worker = WorkerThread(cmd)
        self.worker.finished.connect(self.visualization_finished)
        self.worker.start()
    
    def visualization_finished(self, success, output):
        if success:
            self.status_label.setText("Status: Visualization generated successfully")
            selected_file = self.ms_file_combo.currentText()
            base_name = os.path.splitext(selected_file)[0]
            
            # Look for the newly created PNG file
            possible_files = [f for f in os.listdir('.') if f.startswith(base_name) and f.endswith('.png')]
            
            if possible_files:
                self.display_image(possible_files[0])
            else:
                QMessageBox.warning(self, "Warning", "Visualization generated but no output image found.")
        else:
            self.status_label.setText("Status: Error generating visualization")
            QMessageBox.critical(self, "Error", f"Error generating visualization:\n{output}")
        
        self.generate_btn.setEnabled(True)
    
    def display_image(self, image_path):
        pixmap = QPixmap(image_path)
        if not pixmap.isNull():
            self.image_label.setPixmap(pixmap)
            self.image_label.resize(pixmap.size())
            self.image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        else:
            self.image_label.setText(f"Cannot load image: {image_path}")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = EdacGUI()
    window.show()
    sys.exit(app.exec()) 