try:
    import pytcpb as tc
except ImportError:
    raise ImportError("Failed to import pytcpb")
import sys
from ase.io import read
from ase import units
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--tcfile", type=str, default="s0.inp", help="TeraChem input file")
parser.add_argument("--xyzfile", type=str, help="xyz file to read in")
parser.add_argument("--port", type=int, default=8080, help="TeraChem server port")

args = parser.parse_args()


# Set information about the server
host = "localhost"

# Get Port
port = args.port

# Example TC Input
tcfile = args.tcfile

# Set global treatment (for how TeraChem will handle wavefunction initial guess)
# 0 means continue and use casguess, 1 is keep global variables the same, but recalc casguess, 2 means reinitialize everything
globaltreatment = {"Cont": 0, "Cont_Reset": 1, "Reinit": 2}

# Information about initial QM region
structures = read(args.xyzfile, index=":")
qmattypes = structures[0].get_chemical_symbols()

# Attempts to connect to the TeraChem server
print(f"Attempting to connect to TeraChem server using host {host} and {port}.")
status = tc.connect(host, port)
if status == 0:
    print("Connected to TC server.")
elif status == 1:
    raise ValueError("Connection to TC server failed.")
elif status == 2:
    raise ValueError(
        "Connection to TC server succeeded, but the server is not available."
    )
else:
    raise ValueError("Status on tc.connect function is not recognized!")


# Setup TeraChem
status = tc.setup(str(tcfile), qmattypes)
if status == 0:
    print("TC setup completed with success.")
elif status == 1:
    raise ValueError(
        "No options read from TC input file or mismatch in the input options!"
    )
elif status == 2:
    raise ValueError("Failed to setup TC.")
else:
    raise ValueError("Status on tc_setup function is not recognized!")

for i, structure in enumerate(structures):
    qmcoords = structure.get_positions().flatten() / units.Bohr
    qmcoords = qmcoords.tolist()

    # Compute energy and gradient
    time.sleep(0.010)  # TCPB needs a small delay between calls
    if i == 0:
        totenergy, qmgrad, mmgrad, status = tc.compute_energy_gradient(
            qmattypes, qmcoords, globaltreatment=globaltreatment["Cont_Reset"]
        )
        print("Starting new positions file")
    else:
        print(f"Continuing with job {i}")
        totenergy, qmgrad, mmgrad, status = tc.compute_energy_gradient(
            qmattypes, qmcoords, globaltreatment=globaltreatment["Cont"]
        )

    print(f"Status: {status}")
    if status == 0:
        print("Successfully computed energy and gradients")
    elif status == 1:
        raise ValueError("Mismatch in the variables passed to compute_energy_gradient")
    elif status == 2:
        raise ValueError("Error in compute_energy_gradient.")
    else:
        raise ValueError(
            "Status on compute_energy_gradient function is not recognized!"
        )
