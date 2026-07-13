import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import os
import numpy as np
import matplotlib.pyplot as plt

class newtonspecopener():
    def __init__(self, root, opennotebook):
        self.root = root
        self.opennotebook = opennotebook
        self.buildopenframe()
    
    def buildopenframe(self):
        self.openframe = tk.Frame(self.root, bg = "white")

