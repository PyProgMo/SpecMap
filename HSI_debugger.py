class debug_main:
    def __init__(self, main): #debug Main application
        self.main = main

class debug_Hyperspectral:
    def __init__(self, hyperspectral): #debug Hyperspectral
        self.hyperspectral = hyperspectral

class debug_Clara_Image:
    def __init__(self, clara_image): #debug Clara Image
        self.clara_image = clara_image

class debug_Export:
    def __init__(self, export): #debug Export
        self.export = export

class debug_Newton_Spectrum:
    def __init__(self, newton_spectrum): #debug Newton Spectrum
        self.newton_spectrum = newton_spectrum

class debug_Settings:
    def __init__(self, settings): #debug Settings
        self.settings = settings

class debug_HSI: # debug Hyperspectral frame
    def __init__(self, hsi):
        self.hsi = hsi

class main_Debugger:
    def __init__(self):
        pass

    def build_GUI(self, maingui):
        self.maingui = maingui # maingui must be tkinter root element to place GUI
        self.debug_windows = {}
        self.build_GUI()

    def construct_debug_windows(self):
        self.debug_windows = {
            "main": debug_main(self.main),
            "hyperspectral": debug_Hyperspectral(self.main.hyperspectral),
            "clara_image": debug_Clara_Image(self.main.clara_image),
            "export": debug_Export(self.main.export),
            "newton_spectrum": debug_Newton_Spectrum(self.main.newton_spectrum),
            "settings": debug_Settings(self.main.settings),
            "hsi": debug_HSI(self.main.hsi),
        }
    
    def build_GUI(self):
        self.construct_debug_windows()
        # Build the GUI elements for each debug window
        for name, debug_window in self.debug_windows.items():
            debug_window.build_GUI(self.maingui)
