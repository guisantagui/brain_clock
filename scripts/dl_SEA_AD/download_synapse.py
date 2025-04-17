################################################################################
# Download files or directories from synapse.org                               #
################################################################################

import synapseclient
from synapseutils import syncFromSynapse
import argparse
import os
from dotenv import load_dotenv

# Terminal argument parser
parser = argparse.ArgumentParser(description='Download files from synapse.org given a synID.')

parser.add_argument('--synID', 
	                  type=str,
                    help='Synapse ID with the file.')

parser.add_argument('--outDir', 
	                  type=str,
                    help='Directory where files will be stored.', 
                    default=os.getcwd())

args = parser.parse_args()

# Directory stuff
################################################################################
synID = args.synID
outDir = args.outDir

if not os.path.isdir(outDir):
    os.makedirs(outDir, exist_ok=True)

# Functions
################################################################################
def download_synapse_entity(syn_id, download_path):
    # Get entity metadata
    entity = syn.get(syn_id, downloadFile=False)
    
    # If entity is a file, download directly
    if entity.concreteType == "org.sagebionetworks.repo.model.FileEntity":
        file_entity = syn.get(
            syn_id,
            downloadLocation=download_path,
            ifcollision = 'keep.local'
        )
        print(f"Downloaded file to: {file_entity.path}")

    # If entity is a folder or a project, use syncFromSynapse
    elif entity.concreteType in [
        "org.sagebionetworks.repo.model.Folder",
        "org.sagebionetworks.repo.model.Project"
    ]:
        syncFromSynapse(
            syn,
            syn_id,
            path=download_path,
            ifcollision = 'keep.local'
        )
        print(f"Downloaded folder contents to: {download_path}")
    else:
        print("Unsupported Synapse entity type!")

# Authenticate
################################################################################
load_dotenv()
token = os.getenv("ACCESS_TOKEN")

syn = synapseclient.Synapse()
syn.login(authToken = token)

# Download file
################################################################################
download_synapse_entity(synID, outDir)