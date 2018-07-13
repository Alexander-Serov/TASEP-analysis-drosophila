

import os
import shutil


def reinit_folder(folder):

    if os.path.isdir(folder):
        try:
            shutil.rmtree(folder)
        except Exception as e:
            print(e)

    # Recreate the folder
    try:
        os.makedirs(folder)
    except Exception as e:
        print(e)
