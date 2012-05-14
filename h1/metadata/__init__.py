"""
Extract metadata from files (e.g. models).

Authors: Bernhard Seiwald <bernhard.seiwald@gmail.com>
         Dave Pretty <david.pretty@anu.edu.au>



"""
import mimetypes
from h1.metadata import filecheck

supported_mimetypes = {
    'application/x-netcdf':filecheck.netcdf,
    }

def read_metadata(filename):
    """
    ...
    """
    # Raise IOError if file doesn't exist.
    # We do this as a separate step as guess_mimemtype
    # returns None (less useful) if the file doesn't exist.
    with open(filename) as f: pass

    # determine mimetype
    file_mimetype = mimetypes.guess_type(filename)[0]
        
    if not file_mimetype in supported_mimetypes.keys():
        raise NotImplementedError("File type not supported.")
    
    return supported_mimetypes[file_mimetype](filename)
