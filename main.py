# import numpy as np
# import pyvista as pv
# import panel as pn
# import param
# from scipy.constants import mu_0

# # Inicjalizacja rozszerzenia VTK dla Panel
# pn.extension('vtk', sizing_mode="stretch_both")

# class EMSimulator(param.Parameterized):
#     # --- Parametry Prądu ---
#     I0 = param.Number(5.0, bounds=(-20, 20), label="Natężenie I0 [A]")
#     L = param.Number(2.0, bounds=(0.5, 5.0), label="Długość L [m]")
#     r_wire = param.Number(0.05, bounds=(0.01, 0.2), label="Promień r [m]")
#     current_type = param.ObjectSelector(default="Stały (DC)", objects=["Stały (DC)", "Zmienny (AC)"])
#     freq = param.Number(50.0, bounds=(0.1, 1000.0), label="Częstotliwość f [Hz]")
#     time_phase = param.Number(0.0, bounds=(0, 2*np.pi), label="Faza (ωt) [rad]")

#     # --- Wizualizacja ---
#     show_B = param.Boolean(False, label="Pokaż pole magnetyczne (B)")
#     show_E = param.Boolean(False, label="Pokaż pole elektryczne (E)")
#     show_wire = param.Boolean(True, label="Pokaż przewodnik")
    
#     # --- Kamera ---
#     cam_x = param.Number(5.0, bounds=(-10, 10), label="Kamera X")
#     cam_y = param.Number(5.0, bounds=(-10, 10), label="Kamera Y")
#     cam_z = param.Number(5.0, bounds=(-10, 10), label="Kamera Z")
#     zoom = param.Number(1.0, bounds=(0.1, 5.0), label="Przybliżenie (Zoom)")

#     def __init__(self, **params):
#         super().__init__(**params)
#         # Tworzymy plotter, ale nie renderujemy nic na starcie
#         self.plotter = pv.Plotter(border=False)
#         self.plotter.background_color = "#1e1e1e"
#         self.vtk_pane = pn.pane.VTK(self.plotter.ren_win, width=900, height=900, sizing_mode='stretch_both')

#     def calculate_fields(self):
#         """Obliczenia fizyczne: Biot-Savart i Faraday"""
#         res = 25
#         x_s = np.linspace(-self.L, self.L, res)
#         y_s = np.linspace(-self.L, self.L, res)
#         z_s = np.linspace(-self.L, self.L, res)
#         x, y, z = np.meshgrid(x_s, y_s, z_s, indexing='ij')
        
#         pts = np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)
#         rho = np.sqrt(pts[:,0]**2 + pts[:,1]**2) + 1e-6 # Epsilon
#         phi = np.arctan2(pts[:,1], pts[:,0])
        
#         # Prąd chwilowy
#         I_t = self.I0 if self.current_type == "Stały (DC)" else self.I0 * np.sin(self.time_phase)
        
#         # Pole B (Analityczne rozwiązanie Biota-Savarta dla odcinka)
#         z1, z2 = -self.L/2, self.L/2
#         d1 = np.sqrt(rho**2 + (pts[:,2] - z1)**2)
#         d2 = np.sqrt(rho**2 + (pts[:,2] - z2)**2)
#         B_mag = (mu_0 * I_t / (4 * np.pi * rho)) * ((pts[:,2] - z1)/d1 - (pts[:,2] - z2)/d2)
        
#         B_vec = np.zeros_like(pts)
#         B_vec[:, 0] = -B_mag * np.sin(phi)
#         B_vec[:, 1] = B_mag * np.cos(phi)
        
#         # Pole E (Indukcja dla AC)
#         E_vec = np.zeros_like(pts)
#         if self.current_type == "Zmienny (AC)":
#             omega = 2 * np.pi * self.freq
#             didt = self.I0 * omega * np.cos(self.time_phase)
#             # Czynnik geometryczny dla E_phi
#             f_z = (pts[:,2]+self.L/2)/np.sqrt(rho**2+(pts[:,2]+self.L/2)**2) - \
#                   (pts[:,2]-self.L/2)/np.sqrt(rho**2+(pts[:,2]-self.L/2)**2)
#             E_mag = -(mu_0 * didt / (4 * np.pi)) * f_z
#             E_vec[:, 0] = -E_mag * np.sin(phi)
#             E_vec[:, 1] = E_mag * np.cos(phi)
            
#         return pts, B_vec, E_vec

#     @param.depends('I0', 'L', 'r_wire', 'current_type', 'freq', 'time_phase', 
#                    'show_B', 'show_E', 'show_wire', 'cam_x', 'cam_y', 'cam_z', 'zoom', watch=True)
#     def update_view(self):
#         self.plotter.clear()
        
#         # Aktualizacja kamery
#         self.plotter.camera_position = [(self.cam_x, self.cam_y, self.cam_z), (0,0,0), (0,0,1)]
#         self.plotter.camera.zoom(self.zoom)

#         if self.show_wire:
#             wire = pv.Cylinder(center=(0,0,0), direction=(0,0,1), radius=self.r_wire, height=self.L)
#             self.plotter.add_mesh(wire, color="orange", opacity=0.3)

#         pts, B, E = self.calculate_fields()
#         grid = pv.PolyData(pts)

#         if self.show_B:
#             grid["B"] = B
#             arrows = grid.glyph(orient="B", scale="B", factor=0.8, geom=pv.Arrow())
#             self.plotter.add_mesh(arrows, cmap="viridis", label="Pole B")

#         if self.show_E and self.current_type == "Zmienny (AC)":
#             grid["E"] = E
#             arrows_e = grid.glyph(orient="E", scale="E", factor=0.5, geom=pv.Arrow())
#             self.plotter.add_mesh(arrows_e, color="cyan", label="Pole E")

#         self.vtk_pane.param.trigger('object') # Wymuszenie odświeżenia widżetu

#     def action_reset_view(self, event):
#         self.cam_x, self.cam_y, self.cam_z, self.zoom = 5.0, 5.0, 5.0, 1.0

# # --- UI Setup ---
# sim = EMSimulator()

# reset_btn = pn.widgets.Button(name="Zresetuj Widok", button_type="warning")
# reset_btn.on_click(sim.action_reset_view)

# update_btn = pn.widgets.Button(name="Uruchom / Zaktualizuj Symulację", button_type="success")
# update_btn.on_click(lambda e: sim.update_view())

# control_panel = pn.Column(
#     "# Panel Sterowania",
#     pn.Tabs(
#         ("Parametry", pn.Column(
#             "### Prąd i Geometria",
#             sim.param.I0, sim.param.L, sim.param.r_wire,
#             sim.param.current_type, sim.param.freq, sim.param.time_phase
#         )),
#         ("Wizualizacja", pn.Column(
#             "### Widoczność",
#             sim.param.show_B, sim.param.show_E, sim.param.show_wire,
#             update_btn
#         )),
#         ("Kamera", pn.Column(
#             "### Pozycja Kamery",
#             sim.param.cam_x, sim.param.cam_y, sim.param.cam_z, sim.param.zoom,
#             reset_btn
#         ))
#     ),
#     width=450,
#     styles={'background': '#2b2b2b', 'padding': '15px'}
# )

# # Layout 16:9
# layout = pn.Row(
#     sim.vtk_pane,
#     control_panel,
#     sizing_mode='stretch_both',
#     styles={'background': '#1e1e1e'}
# )

# layout.show(title="Fizyka Pola EM - Symulacja 3D")

import sys
import numpy as np
import pyvista as pv
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QSlider, QLabel, QCheckBox, 
                             QSplitter, QFrame, QGroupBox, QScrollArea, 
                             QDoubleSpinBox)
from PyQt5.QtCore import Qt, QTimer, pyqtSignal
from pyvistaqt import QtInteractor

# --- Stałe Fizyczne ---
MU_0 = 4 * np.pi * 1e-7  # Przenikalność magnetyczna próżni (T*m/A)
EPS_0 = 8.854e-12        # Przenikalność elektryczna próżni (F/m)

class SimulationState:
    """
    Klasa przechowująca stan całej symulacji.
    Oddziela dane (Model) od interfejsu (View).
    """
    def __init__(self):
        # Parametry Przewodnika 1
        self.w1_curr = 50.0      # Natężenie [A]
        self.w1_pos_x = -1.50    # Pozycja X [m]
        self.w1_len = 1.40       # Długość [m]
        
        # Parametry Przewodnika 2
        self.w2_curr = 50.0
        self.w2_pos_x = 1.50
        self.w2_len = 6.00
        
        # Parametry Pola
        self.charge_density = 0.0  # Gęstość ładunku [C/m]
        self.show_b = True         # Pokaż pole B
        self.show_e = False        # Pokaż pole E
        
        # Animacja
        self.anim_speed = 0.10
        self.is_playing = False
        self.time = 0.0            # Czas symulacji (faza)

class FloatSlider(QWidget):
    """
    Widget łączący suwak (Slider) z polem edycji (SpinBox).
    Pozwala na zmianę wartości myszką lub wpisanie jej z klawiatury.
    """
    valueChanged = pyqtSignal(float)

    def __init__(self, min_val, max_val, init_val, label_text, resolution=100, decimals=2):
        super().__init__()
        self.resolution = resolution
        self.val = init_val
        
        # Główny układ
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 5, 0, 5)
        layout.setSpacing(2)
        
        # Górny pasek: Etykieta nazwy i Pole edycji (SpinBox)
        top_layout = QHBoxLayout()
        self.name_label = QLabel(label_text)
        
        # Konfiguracja pola do wpisywania liczb
        self.spinbox = QDoubleSpinBox()
        self.spinbox.setRange(min_val, max_val)       # Zakres taki sam jak suwaka
        self.spinbox.setValue(init_val)               # Wartość początkowa
        self.spinbox.setDecimals(decimals)            # Ilość miejsc po przecinku
        self.spinbox.setSingleStep(10 / resolution)   # Krok przy kliknięciu strzałek
        self.spinbox.setKeyboardTracking(False)       # Aktualizuj dopiero po Enter/FocusOut (opcjonalne)
        self.spinbox.setFixedWidth(80)                # Stała szerokość dla estetyki

        # Połączenie sygnału zmiany wartości w polu tekstowym
        self.spinbox.valueChanged.connect(self._on_spinbox_change)
        
        top_layout.addWidget(self.name_label)
        top_layout.addStretch()
        top_layout.addWidget(self.spinbox)
        layout.addLayout(top_layout)
        
        # Konfiguracja Suwaka
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMinimum(int(min_val * resolution))
        self.slider.setMaximum(int(max_val * resolution))
        self.slider.setValue(int(init_val * resolution))
        self.slider.valueChanged.connect(self._on_slider_change)
        layout.addWidget(self.slider)
        
        self.setLayout(layout)

    def _on_slider_change(self, int_val):
        """Wywoływane, gdy ruszasz suwakiem."""
        float_val = int_val / self.resolution
        self.val = float_val
        
        # Aktualizujemy SpinBox, ale BLOKUJEMY jego sygnały, 
        # żeby nie wywołał zwrotnie funkcji _on_spinbox_change
        self.spinbox.blockSignals(True)
        self.spinbox.setValue(float_val)
        self.spinbox.blockSignals(False)
        
        self.valueChanged.emit(float_val)

    def _on_spinbox_change(self, float_val):
        """Wywoływane, gdy wpisujesz wartość ręcznie."""
        self.val = float_val
        
        # Aktualizujemy Suwak, blokując jego sygnały
        self.slider.blockSignals(True)
        self.slider.setValue(int(float_val * self.resolution))
        self.slider.blockSignals(False)
        
        self.valueChanged.emit(float_val)

    def set_value(self, val):
        """Zewnętrzne ustawienie wartości."""
        self.slider.setValue(int(val * self.resolution))
        self.spinbox.setValue(val)
        self.val = val

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Symulacja Pól Elektromagnetycznych - Native PyVista & Qt")
        self.resize(1400, 900) # Format 16:9 przybliżony

        # Inicjalizacja Stanu
        self.state = SimulationState()

        # Główny Widget (Splitter)
        self.splitter = QSplitter(Qt.Horizontal)
        self.setCentralWidget(self.splitter)

        # 1. Lewy Panel (ScrollArea dla kontrolek)
        self.scroll = QScrollArea()
        self.control_panel = QWidget()
        self.control_layout = QVBoxLayout()
        self.control_panel.setLayout(self.control_layout)
        self.scroll.setWidget(self.control_panel)
        self.scroll.setWidgetResizable(True)
        self.scroll.setMinimumWidth(300)
        self.splitter.addWidget(self.scroll)

        # 2. Prawy Panel (Wizualizacja PyVista)
        self.frame = QFrame()
        self.plotter = QtInteractor(self.frame)
        vlayout = QVBoxLayout()
        vlayout.setContentsMargins(0, 0, 0, 0)
        vlayout.addWidget(self.plotter.interactor)
        self.frame.setLayout(vlayout)
        self.splitter.addWidget(self.frame)
        
        # --- POPRAWKA: Jawne ustawienie rozmiarów paneli [lewy, prawy] ---
        self.splitter.setSizes([350, 1050])

        # Budowa Interfejsu
        self.setup_ui()
        
        # Inicjalizacja Siatki i Sceny 3D
        self.init_3d_scene()
        
        # Timer Animacji
        self.timer = QTimer()
        self.timer.timeout.connect(self.animation_loop)

        # Pierwsze Rysowanie
        self.update_simulation()

    def setup_ui(self):
        """Tworzy grupy kontrolek zgodnie z wymaganiami użytkownika."""
        
        # --- Grupa: Przewodnik 1 ---
        g1 = QGroupBox("Przewodnik 1")
        l1 = QVBoxLayout()
        
        self.s_w1_curr = FloatSlider(-100, 100, self.state.w1_curr, "W1 curr (A)")
        self.s_w1_curr.valueChanged.connect(lambda v: self.update_param('w1_curr', v))
        l1.addWidget(self.s_w1_curr)
        
        self.s_w1_pos = FloatSlider(-5, 5, self.state.w1_pos_x, "W1 pos x (m)")
        self.s_w1_pos.valueChanged.connect(lambda v: self.update_param('w1_pos_x', v))
        l1.addWidget(self.s_w1_pos)
        
        self.s_w1_len = FloatSlider(0.1, 10, self.state.w1_len, "W1 len (m)")
        self.s_w1_len.valueChanged.connect(lambda v: self.update_param('w1_len', v))
        l1.addWidget(self.s_w1_len)
        
        g1.setLayout(l1)
        self.control_layout.addWidget(g1)

        # --- Grupa: Przewodnik 2 ---
        g2 = QGroupBox("Przewodnik 2")
        l2 = QVBoxLayout()
        
        self.s_w2_curr = FloatSlider(-100, 100, self.state.w2_curr, "W2 curr (A)")
        self.s_w2_curr.valueChanged.connect(lambda v: self.update_param('w2_curr', v))
        l2.addWidget(self.s_w2_curr)
        
        self.s_w2_pos = FloatSlider(-5, 5, self.state.w2_pos_x, "W2 pos x (m)")
        self.s_w2_pos.valueChanged.connect(lambda v: self.update_param('w2_pos_x', v))
        l2.addWidget(self.s_w2_pos)
        
        self.s_w2_len = FloatSlider(0.1, 10, self.state.w2_len, "W2 len (m)")
        self.s_w2_len.valueChanged.connect(lambda v: self.update_param('w2_len', v))
        l2.addWidget(self.s_w2_len)
        
        g2.setLayout(l2)
        self.control_layout.addWidget(g2)

        # --- Grupa: Parametry Pola ---
        gf = QGroupBox("Parametry Pola")
        lf = QVBoxLayout()
        
        # Ustawiamy decimals=9, żeby widzieć wartości rzędu nano-coulombów
        self.s_charge = FloatSlider(-1e-6, 1e-6, self.state.charge_density, "Charge density (C/m)", resolution=10000000, decimals=9)
        self.s_charge.valueChanged.connect(lambda v: self.update_param('charge_density', v))
        lf.addWidget(self.s_charge)
        
        self.cb_show_b = QCheckBox("Show B (Magnetic)")
        self.cb_show_b.setChecked(self.state.show_b)
        self.cb_show_b.toggled.connect(lambda v: self.update_bool('show_b', v))
        lf.addWidget(self.cb_show_b)
        
        self.cb_show_e = QCheckBox("Show E (Electric)")
        self.cb_show_e.setChecked(self.state.show_e)
        self.cb_show_e.toggled.connect(lambda v: self.update_bool('show_e', v))
        lf.addWidget(self.cb_show_e)
        
        gf.setLayout(lf)
        self.control_layout.addWidget(gf)

        # --- Grupa: Animacja ---
        ga = QGroupBox("Animacja")
        la = QVBoxLayout()
        
        self.s_speed = FloatSlider(0.01, 2.0, self.state.anim_speed, "Anim speed")
        self.s_speed.valueChanged.connect(lambda v: self.update_param('anim_speed', v))
        la.addWidget(self.s_speed)
        
        self.cb_play = QCheckBox("Is playing")
        self.cb_play.setChecked(self.state.is_playing)
        self.cb_play.toggled.connect(self.toggle_animation)
        la.addWidget(self.cb_play)
        
        ga.setLayout(la)
        self.control_layout.addWidget(ga)
        
        self.control_layout.addStretch() # Wypycha kontrolki do góry

    def init_3d_scene(self):
        """Konfiguracja sceny PyVista/VTK."""
        self.plotter.add_axes()
        self.plotter.show_grid()
        self.plotter.set_background("#fafafa") # Jasne tło
        
        # Tworzenie siatki obliczeniowej (Grid)
        # Używamy ImageData dla wydajności (siatka strukturalna)
        # Zakres: x[-4, 4], y[-4, 4], z[-4, 4]
        # Rozdzielczość 20x20x20 punktów (8000 punktów) - optymalne dla realtime
        dim = 25
        self.grid = pv.ImageData()
        self.grid.dimensions = np.array([dim, dim, dim])
        self.grid.origin = (-4, -4, -4)
        self.grid.spacing = (8/(dim-1), 8/(dim-1), 8/(dim-1))
        
        # Kontenery na aktory (żeby móc je usuwać/aktualizować)
        self.actor_vectors = None
        self.actor_wire1 = None
        self.actor_wire2 = None

    def update_param(self, name, value):
        """Aktualizuje stan i odświeża widok (jeśli nie ma animacji)."""
        setattr(self.state, name, value)
        if not self.state.is_playing:
            self.update_simulation()

    def update_bool(self, name, value):
        setattr(self.state, name, value)
        self.update_simulation()

    def toggle_animation(self, is_playing):
        self.state.is_playing = is_playing
        if is_playing:
            # Uruchom timer z interwałem 30ms (~33 FPS)
            self.timer.start(30)
        else:
            self.timer.stop()

    def animation_loop(self):
        """Krok animacji wywoływany przez QTimer."""
        # Zwiększamy czas/fazę
        self.state.time += self.state.anim_speed * 0.1
        self.update_simulation()

    def calculate_field_vectors(self):
        """
        Silnik Fizyczny: Oblicza wektory B i E w każdym punkcie siatki.
        Wykorzystuje wektoryzację NumPy dla maksymalnej wydajności.
        """
        points = self.grid.points
        x = points[:, 0]
        y = points[:, 1]
        z = points[:, 2]

        # Inicjalizacja pól wektorowych
        Bx, By, Bz = np.zeros_like(x), np.zeros_like(x), np.zeros_like(x)
        Ex, Ey, Ez = np.zeros_like(x), np.zeros_like(x), np.zeros_like(x)

        # Funkcja pomocnicza do obliczania pola od jednego przewodu
        def add_wire_field(w_pos_x, w_curr, w_len, phase_offset=0):
            # Prąd zmienny w czasie dla efektu animacji (jeśli Is playing)
            # Jeśli animacja wyłączona, time=0, cos(0)=1 -> prąd stały DC
            # Dodajemy przesunięcie fazy dla ciekawszego efektu
            current_t = w_curr * np.cos(self.state.time + phase_offset)
            
            # Geometria: Przewód wzdłuż osi Z, na pozycji X = w_pos_x, Y = 0
            # Zakres Z przewodu: [-w_len/2, w_len/2]
            
            # Wektory odległości od osi przewodu (dla każdego punktu siatki)
            dx = x - w_pos_x
            dy = y
            rho_sq = dx**2 + dy**2 + 1e-9 # unikanie dzielenia przez 0
            rho = np.sqrt(rho_sq)
            
            # --- POLE MAGNETYCZNE B (Skończony przewód) ---
            if self.state.show_b:
                # Kąty widzenia końców przewodu
                z_start = -w_len / 2
                z_end = w_len / 2
                
                sin_theta2 = (z_end - z) / np.sqrt(rho_sq + (z_end - z)**2)
                sin_theta1 = (z_start - z) / np.sqrt(rho_sq + (z_start - z)**2)
                
                # Wartość skalarna B
                # B = (mu0 * I / 4pi * rho) * (sin_theta2 - sin_theta1)
                # Uwaga: sin_theta1 jest ujemny dla punktów "pomiędzy" końcami
                factor = (MU_0 * current_t) / (4 * np.pi * rho_sq) * (sin_theta2 - sin_theta1)
                
                # Wektor B jest styczny do okręgu wokół przewodu (-y, x, 0)
                # B_vec = factor * rho * vector_product
                # Składowe:
                Bx[:] += -factor * dy
                By[:] += factor * dx
                # Bz jest 0 dla prostego przewodu w osi Z

            # --- POLE ELEKTRYCZNE E (Elektrostatyczne od ładunku) ---
            if self.state.show_e and abs(self.state.charge_density) > 1e-12:
                # E od naładowanego pręta (uproszczenie: pole radialne 1/r)
                # E_vec = (lambda / 2pi eps0 rho) * radial_unit_vec
                # radial_unit_vec = [dx/rho, dy/rho, 0]
                
                # Modulujemy ładunek w czasie jeśli animacja włączona (opcjonalne)
                charge_t = self.state.charge_density # * np.cos(self.state.time)
                
                e_factor = charge_t / (2 * np.pi * EPS_0 * rho_sq)
                
                Ex[:] += e_factor * dx
                Ey[:] += e_factor * dy
                # Dla uproszczenia pomijamy składową Z pola E na końcach pręta

        # Dodaj wpływ Przewodnika 1
        add_wire_field(self.state.w1_pos_x, self.state.w1_curr, self.state.w1_len, phase_offset=0)
        
        # Dodaj wpływ Przewodnika 2 (z przesunięciem fazy PI dla efektu przeciwfazy w animacji)
        add_wire_field(self.state.w2_pos_x, self.state.w2_curr, self.state.w2_len, phase_offset=np.pi)

        # Sumowanie wektorów do wizualizacji
        vectors = np.zeros((len(points), 3))
        
        if self.state.show_b and self.state.show_e:
            # Jeśli oba wybrane, pokazujemy sumę (może być mylące fizycznie, ale ciekawe wizualnie)
            # Lepiej byłoby użyć dwóch kolorów, ale Grid obsługuje jeden wektor aktywny.
            # Znormalizujemy do wizualizacji.
            vectors[:, 0] = Bx * 1e5 + Ex * 1e-3 # Skalowanie arbitralne do widoczności
            vectors[:, 1] = By * 1e5 + Ey * 1e-3
            vectors[:, 2] = Bz * 1e5 + Ez * 1e-3
        elif self.state.show_b:
            vectors[:, 0] = Bx
            vectors[:, 1] = By
            vectors[:, 2] = Bz
        elif self.state.show_e:
            vectors[:, 0] = Ex
            vectors[:, 1] = Ey
            vectors[:, 2] = Ez

        return vectors

    def update_simulation(self):
        """Główna pętla odświeżania widoku."""
        
        # 1. Oblicz wektory
        vectors = self.calculate_field_vectors()
        
        # Magnituda do kolorowania
        mag = np.linalg.norm(vectors, axis=1)
        
        # Zapisz do siatki PyVista
        # Skalowanie logarytmiczne lub pierwiastkowe pomaga wizualizować duże rozpiętości
        self.grid['Vectors'] = vectors
        self.grid['Magnitude'] = mag
        
        # 2. Aktualizacja geometrii przewodników (Cylindry)
        # Usuwamy stare, tworzymy nowe (prosta metoda, przy małej liczbie obiektów jest OK)
        if self.actor_wire1: self.plotter.remove_actor(self.actor_wire1)
        if self.actor_wire2: self.plotter.remove_actor(self.actor_wire2)
        
        # Przewodnik 1
        w1 = pv.Cylinder(center=(self.state.w1_pos_x, 0, 0), 
                         direction=(0, 0, 1), 
                         radius=0.05, 
                         height=self.state.w1_len)
        self.actor_wire1 = self.plotter.add_mesh(w1, color='orange', opacity=0.8)
        
        # Przewodnik 2
        w2 = pv.Cylinder(center=(self.state.w2_pos_x, 0, 0), 
                         direction=(0, 0, 1), 
                         radius=0.05, 
                         height=self.state.w2_len)
        self.actor_wire2 = self.plotter.add_mesh(w2, color='cyan', opacity=0.8)

        # 3. Aktualizacja Pola Wektorowego (Strzałki/Glyphs)
        # Używamy techniki "Glyphing" - stawiamy strzałkę w każdym punkcie siatki
        # Skalujemy strzałki tak, by nie były za wielkie ani za małe
        
        arrows = self.grid.glyph(orient='Vectors', scale='Magnitude', factor=0.2, tolerance=0.05)
        
        if self.actor_vectors:
            self.plotter.remove_actor(self.actor_vectors)
        
        if arrows.n_points > 0:
            self.actor_vectors = self.plotter.add_mesh(arrows, cmap="jet", lighting=True)
            
        # Wymuś odświeżenie renderingu
        self.plotter.update()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    
    # Ustawienie stylu Fusion dla nowoczesnego wyglądu
    app.setStyle("Fusion")
    
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())