#!/usr/bin/env python3

"""
Bulk downloader and extractor for Zenodo MD trajectory datasets.

Workflow
--------
1. Read a list of Zenodo record IDs.
2. Download all archives for each record using zenodo_get.
3. Detect normal archives and split archives.
4. Extract the preferred XTC and TPR.
5. Store them in an organized directory tree.
6. Update catalog.csv.
"""

from __future__ import annotations

import csv
import logging
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import zipfile
import fnmatch

from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Iterable
import argparse

from zenodo_get import download


###############################################################################
#                              USER SETTINGS
###############################################################################

#
# Text file containing one Zenodo deposition ID per line.
#
DEFAULT_RECORD_FILE = Path("records.txt")

#
# Final destination.
#
DEFAULT_OUTPUT_ROOT = Path("output")

#
# Catalog.
#
CATALOG_FILE = Path("catalog.csv")

#
# Log file.
#
LOG_FILE = Path("download_extract.log")

#
# Temporary download directory.
#
DEFAULT_TMP_ROOT = Path(
    os.environ.get(
        "ZENODO_TMPDIR",
        "/dev/shm/zenodo_bulk_extract",
    )
)

#
# Remove downloaded archives after successful extraction.
#
DELETE_ARCHIVES = False

#
# Keep archives that failed inspection/extraction.
#
KEEP_FAILED_ARCHIVES = True


###############################################################################
#                                 LOGGING
###############################################################################


def setup_logging() -> logging.Logger:

    logger = logging.getLogger("zenodo_bulk_extract")

    logger.setLevel(logging.INFO)

    formatter = logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(message)s",
        "%Y-%m-%d %H:%M:%S",
    )

    #
    # Console.
    #

    sh = logging.StreamHandler()

    sh.setFormatter(formatter)

    logger.addHandler(sh)

    #
    # Log file.
    #

    fh = logging.FileHandler(LOG_FILE)

    fh.setFormatter(formatter)

    logger.addHandler(fh)

    return logger


logger = setup_logging()


###############################################################################
#                               UTILITIES
###############################################################################


def ensure_directory(path: Path) -> None:

    path.mkdir(
        parents=True,
        exist_ok=True,
    )


def safe_remove(path: Path) -> None:

    try:

        if path.is_file():

            path.unlink()

        elif path.is_dir():

            shutil.rmtree(path)

    except Exception:

        pass


def member_name(member) -> str:
    """
    Return the filename of a TarInfo or ZipInfo member.
    """

    if isinstance(member, zipfile.ZipInfo):

        return member.filename

    return member.name


def read_record_ids(
    filename: Path,
) -> list[str]:
    """
    Read Zenodo record IDs from a text file.

    Blank lines and lines beginning with '#' are ignored.
    """

    records = []

    with filename.open() as fh:

        for line in fh:

            line = line.strip()

            if not line:

                continue

            if line.startswith("#"):

                continue

            records.append(line)

    return records


def timestamp() -> str:

    return datetime.now().isoformat(
        timespec="seconds",
    )


###############################################################################
#                         ARCHIVE IDENTIFICATION
###############################################################################

#
# Archives we process directly.
#
ARCHIVE_SUFFIXES = (
    ".tar",
    ".tar.gz",
    ".tgz",
    ".zip",
)


#
# Split ZIP fragments.
#
SPLIT_ZIP_PATTERN = ".zip."


###############################################################################
#                             DOWNLOAD
###############################################################################


def download_record(
    record_id: str,
    download_dir: Path,
) -> None:
    """
    Download all archive files belonging to a Zenodo record.

    Non-archive files are ignored by zenodo_get.
    """

    ensure_directory(download_dir)

    logger.info("=" * 72)
    logger.info("Downloading record %s", record_id)
    logger.info("=" * 72)

    download(
        record_or_doi=record_id,
        output_dir=download_dir,
        file_glob=[
            "*.tar",
            "*.tar.gz",
            "*.tgz",
            "*.zip",
            "*.zip.*",
        ],
    )


###############################################################################
#                       ARCHIVE DISCOVERY
###############################################################################


def discover_archives(
    directory: Path,
):
    """
    Discover all archives downloaded into a directory.

    Returns
    -------
    list

    Each element is either

        Path

    for a normal archive,

    or

        {
            "type": "split_zip",
            "base": "...",
            "parts": [...]
        }

    for split ZIP archives.
    """

    archives = []

    #
    # Group split archives by basename.
    #

    split_groups = defaultdict(list)

    for path in sorted(directory.iterdir()):

        if not path.is_file():
            continue

        name = path.name.lower()

        #
        # Split ZIP?
        #

        if ".zip." in name:

            base = name.split(".zip.")[0] + ".zip"

            split_groups[base].append(path)

            continue

        #
        # Ordinary archives.
        #

        if (
            name.endswith(".tar.gz")
            or name.endswith(".tgz")
            or name.endswith(".tar")
            or name.endswith(".zip")
        ):

            archives.append(path)

    #
    # Add split archives.
    #

    for base, parts in split_groups.items():

        parts = sorted(parts)

        archives.append(
            {
                "type": "split_zip",
                "base": base,
                "parts": parts,
            }
        )

    return archives


###############################################################################
#                      SPLIT ZIP UTILITIES
###############################################################################


def is_split_archive(
    archive,
) -> bool:

    return isinstance(archive, dict)


def describe_archive(
    archive,
) -> str:
    """
    Human-readable archive name for logging.
    """

    if is_split_archive(archive):

        return archive["base"]

    return archive.name


def log_discovered_archives(
    archives,
) -> None:
    """
    Log the archives that will be processed.
    """

    if not archives:

        logger.warning(
            "No archives discovered."
        )

        return

    logger.info(
        "Discovered %d archive(s).",
        len(archives),
    )

    for archive in archives:

        if is_split_archive(archive):

            logger.info(
                "  %s (%d split parts)",
                archive["base"],
                len(archive["parts"]),
            )

        else:

            logger.info(
                "  %s",
                archive.name,
            )

###############################################################################
#                   SPLIT ZIP RECONSTRUCTION
###############################################################################

def reconstruct_split_zip(
    archive,
    workdir: Path,
) -> Path:
    """
    Reconstruct a split ZIP archive.

    NOTE
    ----
    This implementation is intentionally left incomplete until we
    experimentally determine the correct reconstruction command for
    zip -s archives.

    Returns
    -------
    Path
        Path to the reconstructed ZIP archive.
    """

    raise NotImplementedError(
        "Split ZIP reconstruction not yet implemented."
    )


def prepare_archive(
    archive,
    workdir: Path,
) -> Path:
    """
    Return a normal archive ready for inspection.

    Ordinary archives are returned unchanged.

    Split ZIP archives are reconstructed first.
    """

    if is_split_archive(
        archive,
    ):

        return reconstruct_split_zip(
            archive,
            workdir,
        )

    return archive



###############################################################################
#                           ARCHIVE INSPECTION
###############################################################################


def archive_members(archive: Path):
    """
    Return the archive object and its members.

    The archive object must remain open while members are used.
    """

    suffix = "".join(archive.suffixes).lower()

    if (
        suffix.endswith(".tar.gz")
        or suffix.endswith(".tgz")
        or suffix.endswith(".tar")
    ):

        tf = tarfile.open(archive, "r:*")

        return tf, tf.getmembers()

    if suffix.endswith(".zip"):

        zf = zipfile.ZipFile(archive)

        return zf, zf.infolist()

    raise RuntimeError(
        f"Unsupported archive: {archive}"
    )


###############################################################################
#                          MEMBER SELECTION
###############################################################################

def find_unique_member(
    members,
    pattern: str,
):
    """
    Find exactly one archive member matching a glob pattern.

    Raises
    ------
    RuntimeError
        If zero or multiple matches are found.
    """

    matches = []

    for member in members:

        name = member_name(member)

        if fnmatch.fnmatch(
            Path(name).name,
            pattern,
        ):

            matches.append(member)

    logger.info(
        "Searching for %s : found %d match(es)",
        pattern,
        len(matches),
    )

    for member in matches:

        logger.info(
            "    %s",
            member_name(member),
        )

    if len(matches) == 0:

        raise RuntimeError(
            f"No file matching '{pattern}'"
        )

    if len(matches) > 1:

        raise RuntimeError(
            f"Multiple files match '{pattern}'"
        )

    logger.info(
        "Selected %s",
        member_name(matches[0]),
    )

    return matches[0]

###############################################################################
#                              EXTRACTION
###############################################################################


def extract_member(
    archive_object,
    member,
    destination: Path,
):
    """
    Extract one member from an archive.
    """

    ensure_directory(
        destination.parent,
    )

    if isinstance(
        archive_object,
        tarfile.TarFile,
    ):

        source = archive_object.extractfile(
            member,
        )

    else:

        source = archive_object.open(
            member,
        )

    if source is None:

        raise RuntimeError(
            "Cannot extract member."
        )

    with source:

        with destination.open("wb") as out:

            shutil.copyfileobj(
                source,
                out,
            )


def inspect_and_extract(
    archive: Path,
    destination: Path,
):
    """
    Inspect an archive and extract the preferred
    XTC and TPR files.

    Returns
    -------
    dict
    """

    logger.info(
        "Inspecting %s",
        archive.name,
    )

    archive_object, members = archive_members(
        archive,
    )

    try:

        xtc = find_unique_member(
            members,
            "*_cluster_center_traj.xtc",
        )

        tpr = find_unique_member(
            members,
            "*_md.tpr",
        )

        xtc_name = Path(
            member_name(xtc)
        ).name

        tpr_name = Path(
            member_name(tpr)
        ).name

        system = xtc_name.removesuffix(
            "_cluster_center_traj.xtc"
        )

        system_dir = (
            destination
            / system
        )

        ensure_directory(
            system_dir,
        )

        xtc_out = (
            system_dir
            / f"{system}.xtc"
        )

        tpr_out = (
            system_dir
            / f"{system}.tpr"
        )

        extract_member(
            archive_object,
            xtc,
            xtc_out,
        )

        extract_member(
            archive_object,
            tpr,
            tpr_out,
        )

        return {

            "system": system,

            "xtc": xtc_out,

            "tpr": tpr_out,

            "archive": archive.name,

            "xtc_member": xtc_name,

            "tpr_member": tpr_name,

        }

    finally:

        archive_object.close()

###############################################################################
#                              CATALOG
###############################################################################


def append_catalog(entry: dict) -> None:
    """
    Append one successful extraction to catalog.csv.
    """

    new_file = not CATALOG_FILE.exists()

    with CATALOG_FILE.open(
        "a",
        newline="",
    ) as fh:

        writer = csv.writer(fh)

        if new_file:

            writer.writerow(
                [
                    "system",
                    "archive",
                    "xtc",
                    "tpr",
                    "timestamp",
                ]
            )

        writer.writerow(
            [
                entry["system"],
                entry["archive"],
                entry["xtc_member"],
                entry["tpr_member"],
                timestamp(),
            ]
        )


###############################################################################
#                          COMMAND LINE
###############################################################################


def parse_args():

    parser = argparse.ArgumentParser(
        description="Bulk download and extract MD trajectories from Zenodo."
    )

    parser.add_argument(
        "records",
        nargs="?",
        type=Path,
        default=DEFAULT_RECORD_FILE,
        help="Text file containing Zenodo record IDs (default: records.txt)",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Output directory",
    )

    parser.add_argument(
        "--tmp",
        type=Path,
        default=DEFAULT_TMP_ROOT,
        help="Temporary download directory",
    )

    parser.add_argument(
        "--keep-archives",
        action="store_true",
        help="Keep downloaded archives after successful extraction",
    )

    return parser.parse_args()

###############################################################################
#                                  MAIN
###############################################################################


def main():

    args = parse_args()
    output_root = args.output
    tmp_root = args.tmp

    ensure_directory(output_root)
    ensure_directory(tmp_root)

    logger.info(
        "Temporary directory: %s",
        tmp_root,
    )

    logger.info(
        "Output directory: %s",
        output_root,
    )

    logger.info(
        "Record list: %s",
        args.records,
    )

    output_root = args.output
    tmp_root = args.tmp

    ensure_directory(output_root)
    ensure_directory(tmp_root)

    record_ids = read_record_ids(
        args.records,
    )

    if not record_ids:

        logger.error(
            "No record IDs found in %s",
            args.records,
        )

        return

    logger.info(
        "Found %d record(s).",
        len(record_ids),
    )

    summary = {
        "records": 0,
        "archives": 0,
        "success": 0,
        "failed": 0,
    }

    for record_id in record_ids:

        summary["records"] += 1

        record_dir = tmp_root / str(record_id)

        ensure_directory(
            record_dir,
        )

        try:

            #
            # Download everything for this record.
            #

            download_record(
                record_id,
                record_dir,
            )

            #
            # Find archives.
            #

            archives = discover_archives(
                record_dir,
            )

            log_discovered_archives(
                archives,
            )

            for archive in archives:

                summary["archives"] += 1

                try:

                    archive_path = prepare_archive(
                        archive,
                        record_dir,
                    )

                    result = inspect_and_extract(
                        archive_path,
                        output_root,
                    )

                    append_catalog(
                        result,
                    )

                    summary["success"] += 1

                    logger.info(
                        "Successfully extracted %s",
                        result["system"],
                    )

                    if not args.keep_archives:

                        safe_remove(
                            archive_path,
                        )

                except Exception:

                    logger.exception(
                        "Failed processing %s",
                        describe_archive(
                            archive,
                        ),
                    )

                    summary["failed"] += 1

        except Exception:

            logger.exception(
                "Failed processing record %s",
                record_id,
            )

    logger.info("=" * 72)
    logger.info("Finished")
    logger.info("=" * 72)
    logger.info(
        "Records   : %d",
        summary["records"],
    )
    logger.info(
        "Archives  : %d",
        summary["archives"],
    )
    logger.info(
        "Succeeded : %d",
        summary["success"],
    )
    logger.info(
        "Failed    : %d",
        summary["failed"],
    )


if __name__ == "__main__":

    main()