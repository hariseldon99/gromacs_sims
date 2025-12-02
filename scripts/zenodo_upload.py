from zenodo_client import Creator, Metadata, ensure_zenodo
from pprint import pprint
import zipfile
import tempfile
import os

homedir = os.environ.get('HOME', '.')  # or wherever you want to store cached data
datadir = os.path.join(homedir, 'gitrepos/gromacs_sims/Antimicrobial_Peptides/Pseudomonas_antibiotics/simulation', 'new_simulation_complexes_param')

complexname = 'ndxgyra2007C'

zipdir = os.path.join(datadir, complexname)
# put the resulting zip file in the system temp directorydir = tempfile.gettempdir()
zipdir = tempfile.gettempdir()
zip_path = os.path.join(zipdir, f'{complexname}.zip')
if os.path.exists(zip_path):
    os.remove(zip_path)

with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
    for root, dirs, files in os.walk(zipdir):
        for fname in files:
            if fname.endswith('.trr'):
                continue
            full_path = os.path.join(root, fname)
            arcname = os.path.relpath(full_path, zipdir)
            zf.write(full_path, arcname)

# Define the metadata that will be used on initial upload
metadata = Metadata(
    title="MD Simulation of DNA gyrase subunit A (gyrA) protein of P. aeruginosa, with mutations recorded in 2007, in complex with Nalidixic Acid (NDX) ligand molecule",
    upload_type='dataset',
    description="""Mutated sequences of DNA gyrase subunit A (gyrA) protein (mutations recorded in 2007) associated with P. aeruginosa are modeled and docked with a molecule of Nalidixic Acid. The data from subsequent Molecular Dynamics simulation is presented.

    The data from the Molecular Modelling and Docking Simulations can be found here:

    Roy, A. (2025). Docking of several proteins associated with strains of Pseudomonas aeruginosa with various antibiotic molecules [Data set]. Zenodo. https://doi.org/10.5281/zenodo.15383904""",
    creators=[
        Creator(
            name='Roy, Analabha',
            affiliation='Department of Physics, The University of Burdwan, India',
            orcid='0000-0002-4797-0624',
        ),
        Creator(
            name='Bandopadhyay, Rajib',
            affiliation='Department of Botany, The University of Burdwan, India',
            orcid='0000-0002-8318-5631',
        ),
    ],
    access_right='embargoed',
    embargo_date='2026-01-01',
    version='1.0',
)

res = ensure_zenodo(
    key=complexname,     # this is a unique key you pick that will be used to store
                        # the numeric deposition ID on your local system's cache
    data=metadata,
    paths=[
        zip_path,
    ],
    sandbox=True,  # remove this when you're ready to upload to real Zenodo
)


pprint(res.json())