#!/usr/bin/env python3

import os
import sys
import argparse
import tarfile
import requests
import logging
from tqdm import tqdm
from requests_toolbelt.multipart.encoder import MultipartEncoder, MultipartEncoderMonitor
from argparse import RawTextHelpFormatter
import math

def create_tarball(directory):
    """Create a .tar.gz archive of `directory` and return its path."""
    # gather all file paths under `directory`
    file_paths = []
    for root, _, files in os.walk(directory):
        for fname in files:
            file_paths.append(os.path.join(root, fname))

    # compute total size for progress bar
    total_size = sum(os.path.getsize(p) for p in file_paths)
    progress = tqdm(total=total_size, unit="B", unit_scale=True, desc="Archiving")
            
    base = os.path.basename(os.path.normpath(directory))
    tarball = f"{base}.tar.gz"
    with tarfile.open(tarball, "w:gz") as tar:
        def _filter(ti):
            if ti.isreg():
                progress.update(ti.size)
            return ti

        for p in file_paths:
            arcname = os.path.join(base, os.path.relpath(p, directory))
            tar.add(p, arcname=arcname, filter=_filter)

    progress.close()
    return tarball

def split_file(file_path, num_parts=20):
    """Split a file into `num_parts` equally sized chunks."""
    total = os.path.getsize(file_path)
    chunk_size = -(-total // num_parts)   # ceiling division
    parts = []
    with open(file_path, "rb") as src:
        for i in range(1, num_parts + 1):
            chunk = src.read(chunk_size)
            if not chunk:
                break
            part_name = f"{file_path}.part{i:02d}"
            with open(part_name, "wb") as dst:
                dst.write(chunk)
            parts.append(part_name)
    return parts


def upload_to_zenodo(token, deposition_id, file_path):
        logging.info(f"Starting upload of {file_path} to deposition {deposition_id}")
        encoder = MultipartEncoder(
            fields={"file": (os.path.basename(file_path), open(file_path, "rb"), "application/gzip")}
        )
        progress = tqdm(total=encoder.len, unit="B", unit_scale=True, desc="Uploading")
        monitor = MultipartEncoderMonitor(encoder, lambda m: progress.update(m.bytes_read - progress.n))
        headers = {
            "Authorization": f"Bearer {token}",
            "Content-Type": monitor.content_type
        }
        url = f"https://zenodo.org/api/deposit/depositions/{deposition_id}/files"
        resp = requests.post(url, data=monitor, headers=headers)
        progress.close()
        resp.raise_for_status()
        logging.info("Upload completed successfully")
        return resp.json()

def main():
    # enable verbose output
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

    
    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,
        description="""Tarball a directory and upload it to Zenodo. Get the deposition ID from the url of the Zenodo web interface as follows:
        1. Create a new publication in Zenodo using the web ui,
        2. Enter some information (e.g., title) and click save,
        3. In the browser, copy the deposition id (e.g., in https://zenodo.org/deposit/12345 , 12345 is the deposition id)
        
        Notes: 
        1. You will need to set the ZENODO_TOKEN environment variable to your Zenodo API token. 
           You can create a token at https://zenodo.org/account/settings/tokens/new.
        2. You may need to update the requests_toolbelt library to the latest version to use MultipartEncoderMonitor. 
           You can do this by running `pip install --upgrade requests_toolbelt`.
        3. The tarball will be created in the current working directory with the same name as the directory being archived.
        """
    )
    parser.add_argument("directory", 
                        help="Directory to archive and then upload or tarball to just upload")
    parser.add_argument("deposition_id", 
                        help="Zenodo deposition ID")
    parser.add_argument(
        "--split",
        type=int,
        nargs="?",
        const=20,
        default=None,
        help="Split the tarball into N parts and upload each part (default: 20)"
    )

    args = parser.parse_args()

    token = os.getenv("ZENODO_TOKEN")
    if not token:
        sys.exit("Error: ZENODO_TOKEN environment variable is not set.")

    if os.path.isdir(args.directory):
        tarball = create_tarball(args.directory)
        print(f"✔ Created archive: {tarball}")
    else:
        tarball = args.directory
        if not os.path.isfile(tarball):
            sys.exit(f"Error: {tarball} is not a valid file or directory.")    
    # After you call create_tarball(...) in main():
    if args.split:
        parts = split_file(tarball, args.split)
        for part in parts:
            print(f"✔ Uploading part: {part}")
            upload_to_zenodo(token, args.deposition_id, part)
    else:    
        result = upload_to_zenodo(token, args.deposition_id, tarball)
        print("✔ Upload successful. Server response:")
        print(result)

if __name__ == "__main__":
    main()
