# add Frame to export PixMatrix
import csv
import tkinter as tk
from tkinter import filedialog
import array_plotclass1 as apc

class Exportframe():
    def __init__(self, Notebook, Nanomap) -> None:
        self.PixMatrix = [[0]]
        self.Notebook = Notebook
        self.Nanomap = Nanomap
        self.gdx = 1
        self.gdy = 1
        self.xaxis = [0]
        self.yaxis = [0]
        self.metadadict = {'Title': 'Pixel Matrix', 'Units': 'Counts', 'xaxis': self.Nanomap.PixAxX, 'yaxis': self.Nanomap.PixAxY}
        self.buildframe()

    def buildframe(self):
        # create a new frame for the image processing
        self.export_frame = tk.Frame(self.Notebook, borderwidth=5, relief="ridge")
        self.export_frame.grid(row=0, column=1, sticky='nsew')
        '''
        self.export_label = tk.Label(self.export_frame, text='Export Pixel Matrix to CSV')
        self.export_label.grid(row=0, column=0)
        self.exportnameentry = tk.Entry(self.export_frame)
        self.exportnameentry.insert(0, 'PixMatrix.csv')
        self.exportnameentry.grid(row=1, column=0)    
        self.export_button = tk.Button(self.export_frame, text='Export Pixel Matrix to .csv format', command=lambda: self.export_to_csv(self.Nanomap.PixMatrix, self.exportname))
        self.export_button.grid(row=2, column=0)
        '''

        # Open a "Save As" dialog to save the file
        self.save_button = tk.Button(self.export_frame, text='Export Pixel Matrix', command=self.save_file)
        self.save_button.grid(row=0, column=0)
    
    def save_file(self):
        file_path = filedialog.asksaveasfilename(
        defaultextension=".csv",  # Default file extension
        filetypes=[("Text files", "*.csv"), ("All files", "*.*")],  # File types to filter
        title="Export PixelMatrix to .csv")  # Title of the dialog
        #self.export_to_csv(self.Nanomap.PixMatrix, file_path)
        if file_path != '':
            self.write_to_csv(self.Nanomap.PixMatrix, file_path)
    
    # export a matrix like the PixMatrix to a csv file
    def export_to_csv(self, matrix, filename='PixMatrix.csv'):
        with open(filename, 'w', newline='') as file:
            writer = csv.writer(file)
            # first row is the header, gdx, gdy, and some further information
            for i in self.metadadict.keys():
                writer.writerow([i, str(self.metadadict[i])])
            # Write the matrix to the CSV file
            writer.writerows(matrix)
    
    def write_to_csv(self, data, filename):
        """
        Writes a 2D array (list of lists) to a .csv file with tab-separated values.

        :param data: List of lists, where each inner list represents a row of data.
        :param filename: The name of the file to write to.
        """
        with open(filename, 'w') as file:
            file.write('Matrix Data\n')
            for i in self.metadadict.keys():
                file.write('{} = {}\n'.format(i, str(self.metadadict[i])))
            for row in data:
                for i in range(len(row)-1):
                    file.write(str(row[i]) + '\t')
                file.write(str(row[-1]) + '\n')

