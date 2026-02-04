from vpython import *
import numpy as np

#scena
scene = canvas(title="Symulacja Pola Magnetycznego (Dynamiczna Siatka)", width=1200, height=700, background=color.black)
scene.forward = vector(-0.8, -0.6, -1)


class MFieldSim():
    def __init__(self):
        self.mi0 = 4 * pi * 1e-7
        self.I = 35
        self.arrows = []
        self.segments = []
        self.segment_coords = []
        self.dist_limit = 5.0
        self.status_text = wtext(text="Aktualny kształt: <b>Okrąg</b>")

    def clear_scene(self):
        for a in self.arrows: a.visible = False
        for s in self.segments: s.visible = False
        self.arrows = []
        self.segments = []
        self.segment_coords = []

    def add_segment(self, start, end):
        self.segment_coords.append((start, end))
        self.segments.append(cylinder(pos=start, axis=end - start, radius=0.08, color=color.orange))

    def calculate_B_at_point(self, P):
        total_B = vector(0, 0, 0)
        min_dist = 999
        epsilon = 0.08

        for start, end in self.segment_coords:
            r1 = P - start
            r2 = P - end
            L_vec = end - start
            L_len = mag(L_vec)
            if L_len == 0: continue

            dist_vec = cross(L_vec, r1)
            d_perp = mag(dist_vec) / L_len
            min_dist = min(min_dist, d_perp)

            if d_perp < epsilon: continue

            u_L = hat(L_vec)
            cos_theta1 = dot(u_L, hat(r1))
            cos_theta2 = dot(-u_L, hat(r2))

            B_dir = hat(cross(L_vec, r1))
            B_mag = (self.mi0 * self.I) / (4 * pi * d_perp) * (cos_theta1 + cos_theta2)
            total_B += B_dir * B_mag

        return total_B, mag(total_B), min_dist

    def update_grid(self, shape_type):
        self.clear_scene()
        self.status_text.text = f"Aktualny kształt: <b>{shape_type}</b>"

        #domyslne ograniczenia (jak nie damy specyficznych)
        z_min, z_max = -3, 3
        x_min, x_max = -self.dist_limit, self.dist_limit

        #ustawienia ograniczajace dla danych ksztaltow
        if shape_type == "Prosty":
            self.add_segment(vector(0, 0, -6), vector(0, 0, 6))
            z_min, z_max = -6, 6

        elif shape_type == "Kształt L":
            self.add_segment(vector(0, 0, -4), vector(0, 0, 2))
            self.add_segment(vector(0, 0, 2), vector(5, 0, 2))
            z_min, z_max = -4, 4
            x_max = 6

        elif shape_type == "Kształt U":
            w, h = 2.5, 4.0
            self.add_segment(vector(w, h, 0), vector(w, -1, 0))
            self.add_segment(vector(w, -1, 0), vector(-w, -1, 0))
            self.add_segment(vector(-w, -1, 0), vector(-w, h, 0))
            z_min, z_max = -3, 3

        elif shape_type == "Okrąg":
            radius = 3.0
            n = 32
            for i in range(n):
                t1, t2 = (2 * pi / n) * i, (2 * pi / n) * (i + 1)
                self.add_segment(vector(radius * cos(t1), radius * sin(t1), 0),
                                 vector(radius * cos(t2), radius * sin(t2), 0))
            z_min, z_max = -3, 3

        #siatka wektorowa xyz
        step = 1.0
        for x in np.arange(x_min, x_max + step, step):
            for y in np.arange(-self.dist_limit, self.dist_limit + step, step):
                for z in np.arange(z_min, z_max + step, step):
                    pos_vec = vector(x, y, z)
                    B_vec, B_mag, d_min = self.calculate_B_at_point(pos_vec)

                    if B_mag > 1e-7:
                        hue = 0.65 * (d_min / self.dist_limit)
                        hue = min(0.65, max(0, hue))

                        self.arrows.append(arrow(
                            pos=pos_vec,
                            axis=hat(B_vec) * 0.6,
                            color=color.hsv_to_rgb(vector(hue, 1, 1)),
                            shaftwidth=0.05,
                            opacity=max(0.3, 1 - (d_min / 8))
                        ))


#nicjalizacja
sim = MFieldSim()
scene.append_to_caption("\n\n Wybierz kształt przewodu: ")


def handle_menu(m):
    sim.update_grid(m.selected)


menu(choices=['Prosty', 'Kształt L', 'Kształt U', 'Okrąg'], index=0, bind=handle_menu)

#punkt startowy od czego zaczynamy
sim.update_grid("Prosty")

while True:
    rate(30)