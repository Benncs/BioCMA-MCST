"""
Script to compress existing HDF5 files using the lossless gzip algorithm.
This allows for faster runtime performance by not using compression during simulations. 
When data is appended, this method enables better chunking and therefore better compression.

Author: CASALE Benjamin 
Date: 10/17/2024
Corresponding version: >0.0.9

Usage:
    This script should be run with the following command:

    python compressor.py <source_file_root> 

    where <source_file_root> is the root path to the existing HDF5 file to be compressed

Example:
    python compressor.py /path/todata 

Dependencies:
    - h5py: For handling HDF5 files
"""

import h5py
import os
import shutil
from typing import List 
import tempfile


def copy_and_compress(h5_source, h5_dest, compression="gzip", compression_opts=9):
    """
    Recursively copies data from h5_source to h5_dest with compression.

    :param h5_source: Source HDF5 file or group
    :param h5_dest: Destination HDF5 file or group
    :param compression: Compression algorithm (default is 'gzip')
    :param compression_opts: Compression level (default is 9 for maximum compression)
    """
    for key in h5_source:
        item = h5_source.get(key, getlink=True)  # Get item with link info

        # Skip copying external links
        if isinstance(item, h5py.ExternalLink):
            print(f"Skipping external link: {key}")
            h5_dest[key] = h5py.ExternalLink(item.filename, item.path)
            continue  # Skip to the next item

        _item = h5_source[key]  # Get the actual item (without link info)

        if isinstance(_item, h5py.Dataset):
            # Copy the dataset
            if _item.shape == ():
                print(f"Copying scalar dataset: {key}")
                h5_dest.create_dataset(key, data=_item[()])
            else:
                print(f"Copying dataset: {key}")
                h5_dest.create_dataset(key, data=_item[()], 
                                       compression=compression, 
                                       compression_opts=compression_opts, 
                                       shuffle=True)

        elif isinstance(_item, h5py.Group):
            # Recursively copy groups
            print(f"Copying group: {key}")
            group = h5_dest.create_group(key)  # Create new group in destination
            copy_and_compress(_item, group, compression=compression, compression_opts=compression_opts)




def scan_hdf5_files(root_dir: str) -> List[str]:
    """
    Scans the given directory and subdirectories for HDF5 files.

    :param root_dir: Directory to scan
    :return: List of absolute paths to HDF5 files
    """
    hdf5_files = []
    for dirpath, _, filenames in os.walk(root_dir):
        for file in filenames:
            if file.endswith('.h5'):
                hdf5_files.append(os.path.join(dirpath, file))
    return hdf5_files

def compress_and_replace(file_path: str, tmp_dir: str):
    """
    Compresses the HDF5 file and replaces the original with the compressed version.

    :param file_path: Path to the original HDF5 file
    :param tmp_dir: Temporary directory to store the compressed file
    """
    tmp_dest = os.path.join(tmp_dir, os.path.basename(file_path))

    with h5py.File(file_path, 'r') as h5_source, h5py.File(tmp_dest, 'w') as h5_dest:
        # Copy and compress all datasets
        copy_and_compress(h5_source, h5_dest)

    # Replace original file with the compressed one
    shutil.move(tmp_dest, file_path)

def parse_args():
    """
    Parses command line arguments.
    """
    import argparse
    parser = argparse.ArgumentParser(description="Compress HDF5 files in a directory using GZIP compression.")
    parser.add_argument('root', type=str, help="Root directory containing HDF5 files to compress.")
    return parser.parse_args()


def main(args):
    # Get all HDF5 files in the directory
    hdf5_files = scan_hdf5_files(args.root)
    
    # Compress and replace each file
    for file in hdf5_files:
        compress_and_replace(file, tempfile.gettempdir())

if __name__ == "__main__":
    args = parse_args()
    main(args)
