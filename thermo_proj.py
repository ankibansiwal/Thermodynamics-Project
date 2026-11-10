import os
import numpy as np
import matplotlib.pyplot as plt

R_BAR = 83.1446
R_J = 8.31446

gas_database = {
    "methane": {"Tc": 190.6, "Pc": 45.99, "omega": 0.012},
    "ethane": {"Tc": 305.3, "Pc": 48.72, "omega": 0.099},
    "propane": {"Tc": 369.8, "Pc": 42.48, "omega": 0.152},
    "n-butane": {"Tc": 425.1, "Pc": 37.96, "omega": 0.200},
    "co2": {"Tc": 304.2, "Pc": 73.83, "omega": 0.224},
    "nitrogen": {"Tc": 126.2, "Pc": 33.98, "omega": 0.038},
    "oxygen": {"Tc": 154.6, "Pc": 50.43, "omega": 0.022},
    "ethylene": {"Tc": 282.3, "Pc": 50.40, "omega": 0.087},
    "ammonia": {"Tc": 405.7, "Pc": 112.80, "omega": 0.253},
    "water": {"Tc": 647.1, "Pc": 220.55, "omega": 0.345}
}

eos_parameters = {
    "VDW": {"sigma": 0.0, "epsilon": 0.0, "Omega": 0.125, "Psi": 27.0 / 64.0},
    "RK": {"sigma": 1.0, "epsilon": 0.0, "Omega": 0.08664, "Psi": 0.42748},
    "SRK": {"sigma": 1.0, "epsilon": 0.0, "Omega": 0.08664, "Psi": 0.42748},
    "PR": {"sigma": 1.0 + np.sqrt(2.0), "epsilon": 1.0 - np.sqrt(2.0), "Omega": 0.07780, "Psi": 0.45724}
}

def get_alpha_terms(model, Tr, omega):
    if model == "VDW":
        return 1.0, 0.0
    elif model == "RK":
        return Tr ** -0.5, -0.5
    elif model == "SRK":
        m = 0.480 + 1.574 * omega - 0.176 * (omega ** 2)
        alpha = (1.0 + m * (1.0 - np.sqrt(Tr))) ** 2
        dln_alpha = -m * np.sqrt(Tr / alpha)
        return alpha, dln_alpha
    elif model == "PR":
        kappa = 0.37464 + 1.54226 * omega - 0.26992 * (omega ** 2)
        alpha = (1.0 + kappa * (1.0 - np.sqrt(Tr))) ** 2
        dln_alpha = -kappa * np.sqrt(Tr / alpha)
        return alpha, dln_alpha

def solve_thermo_state(Tc, Pc, omega, T, P, model="PR"):
    model = model.upper()
    params = eos_parameters[model]
    
    sigma, epsilon = params["sigma"], params["epsilon"]
    Omega, Psi = params["Omega"], params["Psi"]

    Tr = T / Tc
    Pr = P / Pc
    alpha, dln_alpha = get_alpha_terms(model, Tr, omega)

    beta = Omega * (Pr / Tr)
    q = (Psi * alpha) / (Omega * Tr)

    c2 = (epsilon + sigma - 1.0) * beta - 1.0
    c1 = epsilon * sigma * (beta ** 2) - (epsilon + sigma) * beta * (1.0 + beta) + q * beta
    c0 = -(epsilon * sigma * (beta ** 2) * (1.0 + beta) + q * (beta ** 2))

    roots = np.roots([1.0, c2, c1, c0])
    gas_roots = [r.real for r in roots if abs(r.imag) < 1e-7 and r.real > beta]

    if not gas_roots:
        return None

    Z = max(gas_roots)
    V_molar = (Z * R_BAR * T) / P

    if abs(sigma - epsilon) > 1e-7:
        I = (1.0 / (sigma - epsilon)) * np.log((Z + sigma * beta) / (Z + epsilon * beta))
    else:
        I = beta / Z

    HR_RT = (Z - 1.0) + (dln_alpha - 1.0) * q * I
    SR_R = np.log(Z - beta) + dln_alpha * q * I
    ln_phi = (Z - 1.0) - np.log(Z - beta) - q * I

    H_res = HR_RT * R_J * T
    S_res = SR_R * R_J
    phi = np.exp(ln_phi)
    fugacity = phi * P

    return {
        "Model": model,
        "T": T,
        "P": P,
        "Z": Z,
        "V": V_molar,
        "H_R": H_res,
        "S_R": S_res,
        "phi": phi,
        "f": fugacity
    }

def print_results(gas_name, results):
    print(f"\n--- THERMODYNAMIC PROFILE FOR {gas_name.upper()} ---")
    print(f"{'EOS':<6}{'T (K)':<9}{'P (bar)':<10}{'Z':<9}{'V (cm3/mol)':<14}{'H^R (J/mol)':<14}{'S^R (J/mol.K)':<15}{'phi':<9}{'f (bar)':<10}")
    print("-" * 95)
    for r in results:
        print(f"{r['Model']:<6}{r['T']:<9.2f}{r['P']:<10.4f}{r['Z']:<9.4f}{r['V']:<14.2f}{r['H_R']:<14.2f}{r['S_R']:<15.3f}{r['phi']:<9.4f}{r['f']:<10.3f}")
    print("-" * 95)

def generate_plots(gas_name, Tc, Pc, omega, T, model="PR"):
    print(f"\nGenerating {model} thermodynamic plots for {gas_name.upper()} at {T} K...")
    
    pressures = np.linspace(1, 100, 100)
    valid_p, Z_data, phi_data, f_data, H_data, S_data = [], [], [], [], [], []
    
    for p in pressures:
        res = solve_thermo_state(Tc, Pc, omega, T, p, model)
        if res:
            valid_p.append(p)
            Z_data.append(res['Z'])
            phi_data.append(res['phi'])
            f_data.append(res['f'])
            H_data.append(res['H_R'] / 1000.0)
            S_data.append(res['S_R'])
            
    if not valid_p:
        print("No valid gas roots found in this pressure range. Skipping plot.")
        return

    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(f"Thermodynamic Properties of {gas_name.upper()} at {T} K ({model} EOS)", fontsize=16, fontweight='bold')

    axs[0, 0].plot(valid_p, Z_data, 'b-', linewidth=2)
    axs[0, 0].axhline(1, color='gray', linestyle='--', label='Ideal Gas (Z=1)')
    axs[0, 0].set(title="Compressibility Factor vs Pressure", xlabel="Pressure (bar)", ylabel="Z")
    axs[0, 0].grid(True, alpha=0.5)
    axs[0, 0].legend()

    axs[0, 1].plot(valid_p, f_data, 'g-', linewidth=2, label='Real Gas Fugacity')
    axs[0, 1].plot(valid_p, valid_p, 'k--', label='Ideal Gas (f=P)')
    axs[0, 1].set(title="Fugacity vs Pressure", xlabel="Pressure (bar)", ylabel="Fugacity (bar)")
    axs[0, 1].grid(True, alpha=0.5)
    axs[0, 1].legend()

    axs[1, 0].plot(valid_p, H_data, 'r-', linewidth=2)
    axs[1, 0].set(title="Residual Enthalpy vs Pressure", xlabel="Pressure (bar)", ylabel="H^R (kJ/mol)")
    axs[1, 0].grid(True, alpha=0.5)

    axs[1, 1].plot(valid_p, S_data, 'm-', linewidth=2)
    axs[1, 1].set(title="Residual Entropy vs Pressure", xlabel="Pressure (bar)", ylabel="S^R (J/mol*K)")
    axs[1, 1].grid(True, alpha=0.5)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    save_path = os.path.join(script_dir, f"{gas_name}_thermodynamic_plots.png")
    plt.savefig(save_path, dpi=300)
    plt.close()
    
    print(f"Success! Plot saved securely to: {save_path}\n")

def user_menu():
    print("Generalized Cubic Equation of State Calculator")
    print("Available gases:", ", ".join(gas_database.keys()))
    
    while True:
        gas_choice = input("\nEnter gas name (or 'quit' to exit): ").strip().lower()
        if gas_choice == 'quit':
            break
            
        if gas_choice not in gas_database:
            print("Gas not found in database. Try again.")
            continue
            
        try:
            T = float(input("Enter Temperature (K): "))
            P = float(input("Enter Pressure (bar): "))
        except ValueError:
            print("Invalid input. Please use numerical values.")
            continue
            
        eos_choice = input("Enter EOS model (VDW / RK / SRK / PR / ALL) [Default: ALL]: ").strip().upper()
        
        if eos_choice in eos_parameters:
            models_to_run = [eos_choice]
        else:
            models_to_run = ["VDW", "RK", "SRK", "PR"]
            
        props = gas_database[gas_choice]
        
        results = []
        for model in models_to_run:
            res = solve_thermo_state(props["Tc"], props["Pc"], props["omega"], T, P, model)
            if res:
                results.append(res)
        
        if results:
            print_results(gas_choice, results)
        else:
            print("\nError: No valid vapor roots found for these conditions (gas may have condensed).")
            continue
            
        want_plot = input(f"Generate thermodynamic plots for {gas_choice.upper()}? (y/n): ").strip().lower()
        if want_plot == 'y':
            plot_model = "PR" if "PR" in models_to_run else models_to_run[-1]
            generate_plots(gas_choice, props["Tc"], props["Pc"], props["omega"], T, model=plot_model)

if __name__ == "__main__":
    user_menu()
