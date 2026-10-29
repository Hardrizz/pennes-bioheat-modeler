import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class PennesBioheatModeler:
    """
    Biomedical engineering 1D numerical solver for the Pennes Bioheat Equation,
    modeling heat conduction, metabolic heat generation, and blood perfusion in living tissue.
    """

    def __init__(self, length_cm: float = 5.0, nodes: int = 100):
        self.L = length_cm * 1e-2  # Convert cm to meters
        self.N = nodes
        self.dx = self.L / (self.N - 1)
        self.x = np.linspace(0, self.L, self.N)

        # Typical biological tissue properties (Muscle/Liver approximation)
        self.k = 0.52      # Thermal conductivity (W/m·K)
        self.rho = 1050.0  # Tissue density (kg/m^3)
        self.c = 3600.0    # Specific heat capacity of tissue (J/kg·K)
        self.cb = 3840.0   # Specific heat capacity of blood (J/kg·K)
        self.wb = 0.0005   # Blood perfusion rate (kg/m^3·s)
        self.Ta = 37.0     # Arterial blood temperature (°C)
        self.Q_met = 420.0 # Metabolic heat generation (W/m^3)

    def solve_steady_state(self, t_surface_left: float = 37.0, t_surface_right: float = 43.0) -> np.ndarray:
        """
        Solves the steady-state 1D Pennes bioheat equation using finite difference methods.
        Equation: k * (d^2T/dx^2) + wb * cb * (Ta - T) + Q_met = 0
        """
        A = np.zeros((self.N, self.N))
        b = np.zeros(self.N)

        # Left Boundary Condition (Dirichlet)
        A[0, 0] = 1.0
        b[0] = t_surface_left

        # Right Boundary Condition (Dirichlet)
        A[self.N - 1, self.N - 1] = 1.0
        b[self.N - 1] = t_surface_right

        # Interior nodes finite difference setup
        # k * (T_{i+1} - 2*T_i + T_{i-1})/dx^2 - wb*cb*T_i + wb*cb*Ta + Q_met = 0
        coeff = self.k / (self.dx ** 2)
        perfusion_term = self.wb * self.cb

        for i in range(1, self.N - 1):
            A[i, i - 1] = coeff
            A[i, i] = - (2 * coeff + perfusion_term)
            A[i, i + 1] = coeff
            b[i] = - (perfusion_term * self.Ta + self.Q_met)

        # Solve linear system
        temperature_profile = np.linalg.solve(A, b)
        return temperature_profile

    def generate_thermal_plot(self, t_profile: np.ndarray, output_filename: str = 'tissue_temperature_profile.png'):
        """Generates and saves a publication-ready thermal profile plot across the tissue depth."""
        plt.figure(figsize=(9, 5))
        
        distance_mm = self.x * 1000.0  # Convert meters to millimeters
        
        plt.plot(distance_mm, t_profile, color='#e41a1c', linewidth=2.5, label='Tissue Temperature Profile')
        plt.axhline(y=37.0, color='blue', linestyle='--', alpha=0.6, label='Core Body Temp (37°C)')
        plt.axhline(y=42.0, color='orange', linestyle=':', alpha=0.8, label='Thermal Ablation Threshold (42°C)')
        
        plt.title('1D Steady-State Tissue Temperature Distribution (Pennes Bioheat Model)', fontsize=11, fontweight='bold')
        plt.xlabel('Tissue Depth (mm)')
        plt.ylabel('Temperature (°C)')
        plt.grid(True, linestyle=':', alpha=0.6)
        plt.legend(loc='upper left')
        
        plt.tight_layout()
        plt.savefig(output_filename, dpi=300)
        plt.close()
        print(f"[+] Thermal profile plot successfully saved to {output_filename}")


if __name__ == '__main__':
    print("=== PENNES BIOHEAT MODELER EXECUTING ===")
    
    # Initialize model for a 4 cm tissue layer
    # Simulating a hyperthermia treatment where boundary x=L is heated to 45°C
    model = PennesBioheatModeler(length_cm=4.0, nodes=150)
    
    # Solve steady state
    temp_profile = model.solve_steady_state(t_surface_left=37.0, t_surface_right=45.0)

    # Summarize results
    print(f"\n[Thermal Modeling Summary]")
    print(f"  Tissue Thickness     : {model.L * 100:.1f} cm")
    print(f"  Minimum Temperature  : {np.min(temp_profile):.2f} °C")
    print(f"  Maximum Temperature  : {np.max(temp_profile):.2f} °C")
    print(f"  Average Temperature  : {np.mean(temp_profile):.2f} °C")

    # Generate visual report
    model.generate_thermal_plot(temp_profile)
    print("\n[✓] Bioheat simulation complete.")
