import os

class FileDeleter:
    def __init__(self, folder_path, file_type):
        self.folder_path = folder_path
        self.file_type = file_type
        self.count = 0

    def delete_files(self):
        for root, dirs, files in os.walk(self.folder_path):
            for file in files:
                if file.endswith(self.file_type):
                    os.remove(os.path.join(root, file))
                    self.count += 1
                    print(f"Deleted: {os.path.join(root, file)}")

        print(f"Total {self.file_type} files deleted: {self.count}")

#
# folder_path = r"C:\Users\THINKPAD\Desktop\pix2pixRAD\brats_2019\MICCAI_BraTS_2019_Data_Training\HGG"
folder_path = r"C:\Users\THINKPAD\Desktop\pix2pixRAD\brats_2019\MICCAI_BraTS_2019_Data_Training\LGG"

file_type = ".png"
deleter = FileDeleter(folder_path, file_type)
deleter.delete_files()