

import os
import shutil


def reinit_folder(folders):
    """
    Clear folder contents or create the folder if necessary.
    """
    # If string, convert to list
    if isinstance(folders, str):
        folders = [folders]

    for folder in folders:
        if not os.path.isdir(folder):
            os.makedirs(folder)
        else:
            for file in os.listdir(folder):
                file_path = os.path.join(folder, file)
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)

        print(f"Folder '{folder}' cleaned/created successfully!")
