from rdkit import Chem
import os
import requests
import gzip
from io import BytesIO

# Définissez l'URL du fichier SDF que vous souhaitez télécharger
sdf_url = "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_lite_3star.sdf.gz"

# Téléchargez le fichier SDF depuis l'URL
response = requests.get(sdf_url)

# Dossier d'output
output_dir = "output_molecules"
os.makedirs(output_dir, exist_ok=True)

# Vérifiez si le téléchargement a réussi (statut HTTP 200)
if response.status_code == 200:
    # Décompressez le contenu
    decompressed_content = gzip.decompress(response.content)

    # Utilisez RDKit pour lire le contenu décompressé directement depuis la mémoire
    supplier = Chem.SDMolSupplier()
    supplier.SetData(decompressed_content)

    # Iterez à travers les molécules dans le fichier SDF
    for idx, mol in enumerate(supplier):
        if mol is not None:
            # Créez une chaîne contenant les données SDF de la molécule
            sdf_data = Chem.MolToMolBlock(mol)

            # Enregistrez la molécule dans un fichier SDF distinct
            output_filename = f"molecule_{idx + 1}.sdf"
            output_file = os.path.join(output_dir, output_filename)
            
            # Utilisez 'with' pour s'assurer que le writer est correctement fermé
            with Chem.SDWriter(output_file) as writer:
                writer.write(mol)

            print(f"Molécule {idx + 1} enregistrée dans {output_filename}")
else:
    print(f"Échec du téléchargement. Statut HTTP: {response.status_code}")