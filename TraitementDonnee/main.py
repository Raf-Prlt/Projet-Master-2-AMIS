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
output_dir = "output_molecules_v3000"
os.makedirs(output_dir, exist_ok=True)

# Vérifiez si le téléchargement a réussi (statut HTTP 200)
if response.status_code == 200:
    # Décompressez le contenu
    decompressed_content = gzip.decompress(response.content)

    # Utilisez RDKit pour lire le contenu décompressé directement depuis la mémoire
    supplier = Chem.SDMolSupplier()
    supplier.SetData(decompressed_content )

    molecules_with_cycles = [mol for mol in supplier if mol and mol.GetRingInfo().NumRings() > 0]
    # Iterez à travers les molécules dans le fichier SDF
    for idx, mol in enumerate(molecules_with_cycles):
        if mol is not None:
            # Obtenez l'identifiant ChEBI de la molécule en accédant à la propriété spécifique
            chebi_id_with_prefix = mol.GetProp("ChEBI ID")  # Assurez-vous que "ChEBI ID" est le bon nom de champ
            # Retirez le préfixe "CHEBI:" en faisant un split
            chebi_id = chebi_id_with_prefix.split(":")[1] if ":" in chebi_id_with_prefix else chebi_id_with_prefix

            # Créez une chaîne contenant les données SDF de la molécule
            sdf_data = Chem.MolToMolBlock(mol)

            # Enregistrez la molécule dans un fichier SDF distinct avec le nom ChEBI_ID
            output_filename = f"molecule_{chebi_id}.sdf"
            output_file = os.path.join(output_dir, output_filename)
            
            # Utilisez 'with' pour s'assurer que le writer est correctement fermé
            with Chem.SDWriter(output_file) as writer:
                writer.SetForceV3000(True)
                writer.write(mol)

            print(f"Molécule {idx + 1} avec ChEBI ID {chebi_id} enregistrée dans {output_filename}")
else:
    print(f"Échec du téléchargement. Statut HTTP: {response.status_code}")