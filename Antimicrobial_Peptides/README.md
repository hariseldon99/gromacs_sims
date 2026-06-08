
## Antimicrobial Peptides
- Contains subdirectories for various peptide simulations.
- Anti-Microbial Research Association (AMRA), The University of Burdwan 
  - [Professor Rajib Bandopadhyay](https://buruniv.irins.org/profile/196472), Department of Botany 
  - [Dr. Analabha Roy](https://physics.utexas.edu/~daneel), Department of Physics
  - [Dr. Sumit Hira](https://www.sumithira.in/), Department of Zoology
  - [Dr. Rajendra Kr Roy](https://orcid.org/0009-0007-6009-9283), Department of Botany 
  - Dr. Raju Biswas, Department of Botany  
  - Rajdeep Shaw, Department of Botany  
  - [Dr. Rahul Chandra](https://orcid.org/0000-0001-6328-2461), Kazi Nazrul University: Asansol, West Bengal
  - Samrat Daripa, Department of Zoology
  - Argha Nath, Department of Zoology

- Contents:
  - Simulations involving Colicin peptides.
    - [Colpk_phosphatedylethanolamine](Colpk_phosphatedylethanolamine/) ([Raw Data](https://doi.org/10.5281/zenodo.15375411)), 
    - [Colpk_phosphatedylglycerol](Colpk_phosphatedylglycerol/) ([Raw Data](https://doi.org/10.5281/zenodo.15373387)), 
    - [Colpk_pyocyanin](Colpk_pyocyanin) ([Raw Data](https://doi.org/10.5281/zenodo.15335934))
  
  - Simulations involving Colicin PM peptides.
    - [Colpm_phosphatedylethanolamine](Colpm_phosphatedylethanolamine/) ([Raw Data](https://doi.org/10.5281/zenodo.15375478)), 
    - [Colpm_phosphatedylglycerol](Colpm_phosphatedylglycerol/) ([Raw Data](https://doi.org/10.5281/zenodo.15375267)), 
    - [Colpm_pyocyanin](Colpm_pyocyanin) ([Raw Data](https://doi.org/10.5281/zenodo.15354461))
  
  - [DNA_peptide](DNA_peptide/) ([Raw Data](https://doi.org/10.5281/zenodo.15380257)): Simulations involving DNA-peptide interactions.
   
  - Simulations of Ku04AMP01 antimicrobial peptide. ([Raw Data](https://doi.org/10.5281/zenodo.15380036))
    - [Ku04AMP01_linear](Ku04AMP01_linear/) ,
    - [Ku04AMP01_phosphatedylglycerol](Ku04AMP01_phosphatedylglycerol]/)  
    - [Ku04AMP01_phosphatedylethanolamine](Ku04AMP01_phosphatedylethanolamine]/) 
  
  - Simulations of Ku04AMP02 antimicrobial peptide. ([Raw Data](https://doi.org/10.5281/zenodo.19977246))
    - [Ku04AMP02_linear](Ku04AMP02_linear/) ,
    - [Ku04AMP02_phosphatedylglycerol](Ku04AMP02_phosphatedylglycerol]/) 
    - [Ku04AMP02_phosphatedylethanolamine](Ku04AMP02_phosphatedylethanolamine]/) 
  
  - Simulations of Pseudomonas proteins with various antibiotics
    - [Pseudomonas_antibiotics](Pseudomonas_antibiotics) ([Raw Data](https://doi.org/10.5281/zenodo.15383903)) (& links therein)
## Zenodo Community

- Zenodo Community @ [https://zenodo.org/communities/amra/](https://zenodo.org/communities/amra/)

## Logs

- Deletion Note on 20250505
There was a google drive synchroniztion error that led to topology and trajectory data loss. Full restoration was accomplished on this day. Restoration script:
```python
# pip install --upgrade google-api-python-client google-auth-httplib2 google-auth-oauthlib
from __future__ import print_function
import os
import pickle
import datetime
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request

# If modifying these SCOPES, delete the file token.pickle.
SCOPES = ['https://www.googleapis.com/auth/drive']

def authenticate():
    """
    Handles user authentication. If valid credentials are stored in token.pickle, they are used;
    otherwise, the OAuth flow starts using the downloaded credentials.json file.
    """
    creds = None
    # Token file stores the user's credentials.
    if os.path.exists('token.pickle'):
        with open('token.pickle', 'rb') as token:
            creds = pickle.load(token)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            # credentials.json should be your downloaded OAuth 2.0 client credentials.
            flow = InstalledAppFlow.from_client_secrets_file('credentials.json', SCOPES)
            creds = flow.run_local_server(port=0)
        # Save the credentials for the next run.
        with open('token.pickle', 'wb') as token:
            pickle.dump(creds, token)
    return creds

def main():
    # Authenticate and build the Drive API service.
    creds = authenticate()
    service = build('drive', 'v3', credentials=creds)

    # Compute the RFC 3339 timestamp for three months (90 days) ago from now (UTC).
    four_months_ago = datetime.datetime.utcnow() - datetime.timedelta(days=120)
    four_months_ago_str = four_months_ago.isoformat() + 'Z'

    # Build a query to find files in the trash that:
    # - have "trashed" set to true,
    # - whose name contains ".xtc" or ".tpr", and
    # - whose modifiedTime is later than three months ago.
    #
    # Note: Google Drive does not provide a dedicated "deletion date" field,
    # so we use modifiedTime as a proxy. Depending on your usage, this may not perfectly reflect
    # when the file was trashed.
    query = (
        f"trashed = true and "
        #f"(name contains '.xtc' or name contains '.tpr') and "
        f"(name contains '.xtc' or name contains 'colpk') and "
        f"modifiedTime > '{four_months_ago_str}'"
    )

    files_to_restore = []
    page_token = None
    print("Searching for trashed files to restore...\n")
    while True:
        response = service.files().list(
            q=query,
            spaces='drive',
            fields='nextPageToken, files(id, name)',
            pageToken=page_token
        ).execute()
        for file in response.get('files', []):
            print(f"Found file: {file.get('name')} (ID: {file.get('id')})")
            files_to_restore.append(file)
        page_token = response.get('nextPageToken', None)
        if page_token is None:
            break

    if not files_to_restore:
        print("\nNo matching files found in the trash.")
        return

    # Restore each file by updating its 'trashed' property to False.
    print("\nRestoring files...")
    for file in files_to_restore:
        file_id = file.get('id')
        updated_file = service.files().update(
            fileId=file_id,
            body={'trashed': False}
        ).execute()
        print(f"Restored file: {updated_file.get('name')} (ID: {updated_file.get('id')})")

if __name__ == '__main__':
    main()
```
