# MAIN9.PY - SpecMap Application - Hauptmodul

## Überblickd

**main9.py** ist das Hauptmodul der SpecMap-Anwendung und stellt die Benutzeroberfläche (GUI) für die Verarbeitung hyperspektraler Bilddaten (Hyperspectral Imaging - HSI) bereit. Die Anwendung ermöglicht dem Benutzer, verschiedene Arten von spektralen Daten zu laden, zu verarbeiten und zu visualisieren.

**Dateitypen:**
- HSI-Daten (Hyperspectral Images)
- CLARA-Mikroskopibilder
- Newton-Spektren
- TCSPC-Daten (Time-Correlated Single Photon Counting)

---

## Architektur und Komponenten

### Hauptklassen

#### 1. **FileProcessorApp**
Hauptanwendungsklasse, die die komplette GUI verwaltet und alle Datenverarbeitungsprozesse orchestriert.

- **Verantwortung:** GUI-Verwaltung, Datenladen, Verarbeitung, Thread-Management
- **Ort:** Zeile 22 - 734

#### 2. **specfilesorter**
Spezialklasse für Batch-Verarbeitung und intelligente Dateiorganisation von HSI-Daten.

- **Verantwortung:** Dateiscanning, Merging von Folgemessungen, Batch-Processing
- **Ort:** Zeile 745 - 1350

---

## GUI-Struktur und Komponenten

### Notebook-Tabs (Reiter)

Die Anwendung nutzt `ttk.Notebook` zur Organisation in mehreren Tabs:

| Tab-Name | Funktion |
|----------|----------|
| **Load Data** | Dateityp-spezifische Eingabeformulare und Verarbeitungsoptionen |
| **Hyperspectra** | Hauptvisualisierung: Farbkarte und Spektrenplots |
| **Clara Image** | CLARA-Mikroskopibildanzeige |
| **Newton Spectrum** | Newton-Spektrum-Anzeige |
| **TCSPC** | TCSPC-Datenverarbeitung |
| **HSI Plot** | HSI-Exporter für Bildgenerierung |
| **HSI File Sorter** | Batch-Verarbeitung und intelligente Dateiorganisation |

### Load Data Tab - Eingabeformulare

#### A. **Specmap Load Frame** (Zeile 191-227)
Lädt HSI-Daten aus einem Ordner:
- **Ordnerauswahl:** Folder with spectra
- **Dateiname-Filter:** Dateinamen müssen diese Zeichenkette enthalten
- **Dateiformat:** Dateiendung (z.B. `.txt`, `.csv`)
- **Button:** "Load HSI data" trigget `init_spec_loadfiles()`

#### B. **Hintergrund- und Verarbeitungsoptionen** (Zeile 239-334)
```
Background Frame (Links):
├── Multiple Backgrounds (Checkbox)
├── Linear Background (Checkbox)
└── Power Correction (Checkbox)

Cosmic Removal Frame (Rechts):
├── Remove Cosmics (Checkbox)
├── Cosmic Removal Function (Dropdown)
├── Cosmic Threshold (Entry)
├── Cosmic Width (Entry)
├── Laser Spotsize nm (Entry)
├── Calculate 1st Derivative (Checkbox)
├── Calculate 2nd Derivative (Checkbox)
├── Polynom Order Derivative fit (Dropdown)
├── use n points (Entry)
├── norm counts, then derive (Checkbox)
└── norm on I, then derive (Checkbox)
```

#### C. **Clara Load Frame** (Zeile 335-357)
```
├── File Selection (Entry + Browse)
├── Scaling Factor (Dropdown: 20x, 50x, 100x)
└── Button: "Load Clara data"
```

#### D. **Save/Load Frame** (Zeile 359-399)
```
Save Section:
├── Save Path (Entry + Browse)
└── Button: "Save"

Load Section:
├── Load Path (Entry + Browse)
└── Button: "Load"
```

#### E. **Newton Spectrum Frame** (Zeile 401-411)
```
├── File Selection (Entry + Browse)
└── Button: "Load Newton data"
```

#### F. **TCSPC Frame** (Zeile 413-427)
```
├── TCSPC Main Directory (Entry + Browse)
├── Save Directory (Entry)
└── Button: "Load TCSPC data"
```

#### G. **Multiple HSI Processing** (Zeile 269-333)
```
├── Process Multiple HSIs (Checkbox)
├── Save each HSI object (Checkbox)
├── File Main Directory (Entry + Browse)
├── Save Images Directory (Entry + Browse)
└── Save HSI objects Directory (Entry + Browse)
```

### Hyperspectra Tab - Visualisierung

**Komponenten:**
- **Canvas with Scrollbars:** Für größere HSI-Datensätze
- **cmapframe:** Farbkarte der Hyperspektral-Daten
- **specframe:** Spektrenvisualisierung mit Click-to-Spectrum-Funktion
- **Mousewheel-Support:** Scroll-Funktion für einfache Navigation

**Verwandte Methoden:**
- `_on_hyper_canvas_configure()` - Canvas-Resize-Handler
- `_on_hyper_mousewheel()` - Mousewheel-Event-Handler
- `_bind_hyper_mousewheel()` / `_unbind_hyper_mousewheel()` - Event-Binding

---

## Datenverarbeitung und Workflows

### 1. HSI-Daten Laden Workflow

```
Benutzer gibt Ordner/Dateinamen/Dateiformat an
        ↓
spec_loadfiles() aufgerufen
        ↓
Alte Threads & Nanomap-Objekt bereinigen
        ↓
Dateien im angegebenen Ordner scannen
        ↓
XYMap-Objekt erstellen mit gescannten Dateien
        ↓
Cosmic Removal & Derivate berechnen (wenn aktiviert)
        ↓
Exporter initialisieren
        ↓
Farbkarte & Spektren visualisieren
```

**Code-Einstiegspunkt:** `spec_loadfiles()` (Zeile 501-625)

**Wichtige Parameter:**
- `multiple_BG`: Mehrere Hintergrundspektren
- `linearBG`: Lineare Hintergrundsubtraktion
- `removecosmicsBool`: Cosmic-Ray-Entfernung aktivieren
- `cosmicthreshold`: Schwelle für Cosmic-Ray-Erkennung (default: 20)
- `cosmicwidth`: Breite des Cosmic-Ray-Fensters (default: 3)
- `powercorrectionBool`: Energiekorrektur anwenden

### 2. Multiple HSI Processing Workflow

```
Benutzer aktiviert "Process Multiple HSIs"
        ↓
init_spec_loadfiles() (Zeile 593)
        ↓
Für JEDEN Ordner im Hauptverzeichnis:
    ├── Alte Nanomap bereinigen
    ├── spec_loadfiles() für diesen Ordner aufrufen
    ├── buildandPlotIntCmap() mit Speicherung
    └── HSI-Objekt speichern (wenn aktiviert)
```

**Vorteil:** Automatische Batch-Verarbeitung ohne manuelle Intervention

### 3. Datei-Sortierung und Merging Workflow

Siehe Abschnitt **specfilesorter** weiter unten.

---

## FileProcessorApp - Detaillierte Methoden

### Initialisierung
```python
def __init__(self, root, defaults):
    """Initialisiert die gesamte GUI"""
```
- Erstellt Menu-Bar via `createmenue()`
- Erzeugt Notebook-Tabs via `windownotebook()`
- Erstellt "Load Data" Tab-Buttons via `createbuttons()`
- Initialisiert leeresXYMap-Objekt (noch keine Daten geladen)
- Erstellt Exporter für HSI-Bildexporte

### GUI-Aufbau

#### `createmenue()` (Zeile 128-131)
Erstellt die Menüleiste mittels Funktionen aus `deflib`.

#### `windownotebook(notebookentries)` (Zeile 133-143)
```python
Parameter: notebookentries = ['Load Data', 'Hyperspectra', 'Clara Image', ...]
```
Erstellt `ttk.Notebook` und die entsprechenden Tabs.

#### `createbuttons(Notebook)` (Zeile 145-474)
In drei Phasen:
1. **Canvas mit Scrollbars** erstellen (Zeile 145-172)
2. **5 Input-Frames** erstellen und positionieren:
   - open_frame (Specmap)
   - claraloadframe (Clara)
   - saveframe (Save/Load)
   - newtonframe (Newton)
   - tcspcframe (TCSPC)
3. **specfilesorter** initialisieren (Zeile 473)

### Verarbeitungsmethoden

#### `init_spec_loadfiles()` (Zeile 547-592)
**Zweck:** Entry-Point für Multiple HSI Processing oder Einzel-Verarbeitung

**Logik:**
```python
if make_multiple_HSIsbool == True:
    # Loop durch alle Unterordner
    for folder in os.listdir(maindir):
        spec_loadfiles()  # für jeden Ordner
        buildandPlotIntCmap(savetoimage=...)
        saveNanomap(...) wenn aktiviert
else:
    spec_loadfiles()  # einzeln
```

#### `spec_loadfiles()` (Zeile 501-625)
**Zweck:** Hauptdaten-Lade-Methode

**Schritte:**
1. **Thread-Cleanup:** Alle alten Threads stoppen
2. **Nanomap-Cleanup:** Bestehendes Nanomap-Objekt löschen
3. **Frame-Cleanup:** cmapframe und specframe zerstören
4. **Dateien scannen:** os.walk() mit Filter
5. **Frames neu erstellen:** Neue cmapframe und specframe
6. **XYMap erstellen:** Mit gescannten Dateien, Optionen
7. **Verarbeitung:** Power Correction wenn aktiviert
8. **Exporter initialisieren:** Für HSI-Export

**Exceptions verarbeitet:**
- Frame bereits zerstört
- Nanomap-Cleanup-Fehler
- Keine Dateien gefunden

#### `_stop_managed_threads()` (Zeile 657-673)
Thread-Management für sichere Beendigung:
```python
for thread in self.managed_threads:
    if thread.is_alive():
        stop_event.set()
        thread.join(timeout=0.5)
    else:
        threads_to_remove.append(thread)
```

#### `cl_loadfiles()` (Zeile 675-692)
Lädt CLARA-Mikroskopibilder mit Skalierungsfaktor.

#### `newtonloadfiles()` (Zeile 695-698)
Lädt Newton-Spektrendaten.

#### `tcspcloadfiles()` (Zeile 531-543)
Lädt TCSPC-Daten und initialisiert TCSPCprocessor.

### Speichern und Laden

#### `saveNanomap(filename)` (Zeile 703-729)
**Speichert:** Komplettes XYMap-Objekt als Pickle-Datei
```python
✓ Alle HSI-Daten
✓ Verarbeitungsparameter
✓ ROI-Masken
✓ Derivate (wenn berechnet)
```

**Fehlerbehandlung:**
- Leer-Check
- Existenz-Check des Nanomap
- Automatische `.pkl` Extension
- Dateiname-Inkrement bei Duplikaten

#### `loadhsisaved(filename)` (Zeile 732-779)
**Lädt:** Zuvor gespeichertes XYMap-Objekt

**Besonderheiten:**
- Erstellt Nanomap wenn nicht vorhanden
- skip_gui_build=True während Laden
- `update_after_load()` nach dem Laden

### Scroll- und Event-Handler

#### `_on_canvas_configure(event)` (Zeile 475-478)
Passt Canvas-Breite an (min. 500px).

#### `_on_hyper_canvas_configure(event)` (Zeile 480-483)
Passt Hyperspectra-Canvas-Breite an (min. 1200px).

#### `_on_mousewheel(event)` / `_on_hyper_mousewheel(event)`
Behandelt Mousewheel-Scrolling.

#### `_bind_mousewheel(event)` / `_unbind_mousewheel(event)`
Bindet/Entfernt Mousewheel-Events.

### Dateibrowser

#### `browse_save_path()` (Zeile 485-492)
Öffnet Filedialog zum Speichern.

#### `browse_load_path()` (Zeile 494-501)
Öffnet Filedialog zum Laden.

#### `filter_substring(a, b)` (Zeile 654-655)
Hilfsfunktion: Filtert Elemente, die Substring enthalten.

---

## specfilesorter - Batch Verarbeitung

**Zweck:** Intelligente Organisation und Batch-Verarbeitung von HSI-Messdaten über mehrere Tage hinweg.

### Initialisierung

```python
def __init__(self, tkroot, maindir, filename, fileend, savedir, processdir):
    """
    tkroot: Tkinter-Frame für GUI
    maindir: Verzeichnis mit Messdaten
    filename: Dateinamen-Filter
    fileend: Dateiendungs-Filter
    savedir: Zielverzeichnis für sortierte/gemischte Daten
    processdir: Verzeichnis für weitere Verarbeitung
    """
```

### GUI-Aufbau

```
┌─ specfilesorter GUI ─────────────────────┐
│                                           │
│  LEFT PANEL          │   RIGHT PANEL     │
│  ─────────────────── │   ─────────────── │
│  Inputs & Controls   │   Results Treeview│
│  • Main dir browse   │   • Folder        │
│  • Filename filter   │   • Files count   │
│  • File ext filter   │   • Full path     │
│  • Save dir browse   │                   │
│  • Process dir       │   TreeView:       │
│  • Merge checkbox    │   Extended select │
│  • Buttons           │   Scrollable      │
│  • Progress bar      │                   │
│                                           │
└───────────────────────────────────────────┘
```

**Buttons:**
- **Scan:** Scannt maindir nach matchenden Dateien
- **Preview Selected:** Öffnet Ordner im Explorer
- **Process Selected:** Kopiert Dateien & merged
- **Clear:** Löscht Scan-Ergebnisse
- **Cancel:** Bricht laufenden Copy-Prozess ab (nur während Kopieren aktiv)

### Wichtige Methoden

#### `scan_maindir()` (Zeile 1226-1253)
```
os.walk(maindir):
    für jede Datei:
        ├── filename_filter check
        ├── fileend_filter check
        └── append (folder_name, file_count, path)
```
Ergebnis: `self.scan_results = [(name, count, path), ...]`

#### `populate_list()` (Zeile 1255-1258)
Füllt Treeview mit Scan-Ergebnissen.

#### `on_tree_select(event)` (Zeile 1260-1262)
Speichert selektierte Ordner-Pfade in `self.selected_items`.

#### `preview_selected()` (Zeile 1264-1278)
```
für jeden selected_item:
    sys.platform == 'win32':
        os.startfile(path)
    sonst:
        subprocess.run(['xdg-open', path])
```

#### `sort_and_process()` (Zeile 963-1040)
**Logik:**

1. **Duplikat-Handling:** Für jede eindeutige Ordnername nur erste Instanz nehmen
2. **Dateien zuordnen:** Zu jedem Ordner passende Dateien sammeln
3. **Kopieren:** Dateien mit Fortschritts-Update in Zielordner kopieren
4. **Merging:** Wenn `merge_var == 1`, konsekutive Tage zusammenführen

**Thread-Handling:**
```python
if total_files > THRESHOLD (intern):
    # Starten Worker-Thread für Hintergrund-Kopieren
else:
    # Synchrone Kopie (für wenige Dateien)
```

#### `_copy_worker()` (Zeile 1163-1194)
Background-Thread zum Datenkopieren:
- Thread-safe Counter mit Lock
- Stop-Event-Support
- GUI-Updates via `after()`
- Ruft `_merge_consecutive_days()` am Ende

#### `_merge_consecutive_days(savedir)` (Zeile 1041-1161)
**Intelligentes Merging von Folgemessungen**

**Logik:**
```python
1. Regex-Parse: Extrahiert YYYYMMDD-Datum aus Ordnernamen
    Pattern: ^(prefix)(YYYYMMDD)(suffix)$
    
2. Gruppierung: Gruppiert Ordner mit gleichem prefix+suffix
    Schlüssel: (prefix, suffix)
    
3. Sortierung: Sortiert Gruppen nach Datum
    
4. Sequenz-Erkennung: Findet konsekutive Tagesfolgen
    ├── Wenn Datum[i] == Datum[i-1] + 1 Tag:
    │   └── Merge alle Dateien von i nach i-1
    └── Sonst: Neue Sequenz starten
    
5. Konflikt-Handling: Bei Namensduplikaten:
    Suffix mit Datum hinzufügen: image.png → image_20250102.png
```

**Beispiel:**
```
Input-Ordner:
├── Sample_20250101_Run1
├── Sample_20250102_Run1  ← Konsekutiv, mergen
├── Sample_20250103_Run1  ← Konsekutiv, mergen
└── Sample_20250105_Run2  ← Gap (4.1.), neue Sequenz

Output nach Merge:
├── Sample_20250101_Run1  (enthält alle Dateien von 1,2,3)
└── Sample_20250105_Run2
```

#### `_update_progress(value)` (Zeile 1196-1199)
GUI-Thread-safe Progress-Update.

#### `cancel_copy()` (Zeile 1201-1207)
Setzt `_stop_event` zum Abbrechen des Worker-Threads.

#### `_copy_finished()` (Zeile 1209-1233)
Finalisierung nach Kopieren:
- Re-enables alle Buttons
- Zeigt Fortschritt-Nachricht
- Reset der Progress-Bar

#### `clear_list()` (Zeile 1289-1293)
Löscht alle Scan-Ergebnisse und Auswahl.

---

##  Thread-Management

### Thread-Typen

1. **Managed Threads**
   - Liste: `self.managed_threads`
   - Cleanup: `_stop_managed_threads()`
   - Verwendung: Zukünftige Background-Tasks

2. **Copy Worker Thread**
   - `_copy_worker()` (specfilesorter)
   - Daemon-Thread
   - Stop-Event-basierter Abbruch

### Best Practices

```python
# Thread starten
thread = threading.Thread(target=worker_func, daemon=True)
thread.stop_event = threading.Event()  # für Abbruch
thread.start()

# Thread stoppen
thread.stop_event.set()
thread.join(timeout=0.5)  # mit Timeout
```

---

## Abhängigkeiten und Importe

| Modul | Zweck | Zeile |
|-------|-------|-------|
| `tkinter` | GUI-Framework | 1-4 |
| `tkinter.ttk` | Moderne Widgets | 2 |
| `tkinter.filedialog` | Datei-Dialog | 3 |
| `tkinter.messagebox` | Message-Boxen | 3 |
| `matplotlib` | Visualisierung | 5 |
| `PIL` | Bildverarbeitung | 7 |
| `lib9` | XYMap (HSI-Klassifizierung) | 9 |
| `numpy` | Numerik | 10 |
| `deflib1` | Standardwerte & Hilfsfunktionen | 11 |
| `claralib1` | CLARA-Bild-Loader | 12 |
| `export2` | HSI-Export/Bildgenerierung | 13 |
| `newtonspeclib1` | Newton-Spektrum-Loader | 14 |
| `threading` | Multi-Threading | 15 |
| `HSI_debugger` | Debug-Utility | 18 |
| `TCSPClib` | TCSPC-Daten-Loader | 19 |
| `shutil` | Datei-Operationen | 20 |
| `gc` | Garbage Collection | 20 |
| `error_handler` | Fehlerbehandlung & Logging | 21 |
| `datetime` | Datums-/Zeit-Funktionen | 22 |

---

## Fehlerbehandlung

### Error Engine

```python
error_engine = error_handler.ErrorEngine(
    log_file='logs/specmap.log',
    max_bytes=10*1024*1024,  # 10MB
    backup_count=5  # Rollover-Backups
)
```

**Merkmale:**
- Rotating Log-Datei
- Thread-sichere Logging
- Lokalisierte Error-Messages via Tkinter Messagebox
- Zentralisiert in Hauptmodul

### Fehler-Behandlung Strategien

#### Frame-Cleanup
```python
try:
    if hasattr(self, 'cmapframe') and self.cmapframe.winfo_exists():
        self.cmapframe.destroy()
except tk.TclError:
    pass  # Frame bereits zerstört
```

#### Thread-Cleanup
```python
for thread in threading.enumerate():
    if thread.name != 'MainThread':
        stop_event = getattr(thread, 'stop_event', None)
        if stop_event:
            try:
                stop_event.set()
            except Exception as e:
                print(f"Warning: {e}")
```

#### Datei-Operationen
```python
try:
    file_list = os.listdir(folder)
except Exception as e:
    print(f"Error: Could not list directory {folder}: {e}")
    return
```

---

## Startpunkt und Anwendungs-Lebenszyklus

### `__main__` Block (Zeile 1370-1410)

**Initialisierungsschritte:**

1. **Error Engine Setup** (Zeile 1374-1384)
   ```python
   error_engine = error_handler.ErrorEngine(...)
   logger = error_engine.get_logger()
   logger.info("SpecMap application starting")
   ```

2. **Debugger Initialisierung** (Zeile 1387)
   ```python
   debugger = DBG.main_Debugger()
   ```

3. **Standardwerte Laden** (Zeile 1390-1392)
   ```python
   defaults = deflib.initdefaults()
   ```

4. **Tkinter-Fenster Erstellen** (Zeile 1394-1396)
   ```python
   root = tk.Tk()
   root.geometry(f'{width}x{height}')
   ```

5. **App Initialisieren** (Zeile 1398-1401)
   ```python
   app = FileProcessorApp(root, defaults)
   app.error_engine = error_engine  # Weitergabe für Module
   ```

6. **Close-Handler Registrieren** (Zeile 1404)
   ```python
   root.protocol("WM_DELETE_WINDOW", lambda: pressclose(root, app))
   ```

7. **Hauptloop Starten** (Zeile 1408)
   ```python
   root.mainloop()
   ```

### Close Handler

```python
def pressclose(root, app):
    """Sicheres Herunterfahren"""
    # Alle Non-MainThread-Threads als Daemon setzen
    for thread in threading.enumerate():
        if thread is not threading.main_thread():
            thread.daemon = True
    
    root.destroy()
    app.on_closing()
```

---

## Workflow-Diagramme

### HSI Laden - Ablauf

```
┌──────────────────────┐
│ Benutzer:            │
│ Ordner + Filter      │
│ Parameter auswählen  │
└──────────────────────┘
           │
           ↓
┌──────────────────────┐
│ init_spec_loadfiles()│
└──────────────────────┘
           │
      ┌────┴────┐
      │          │
   JA │          │ NEIN
      ↓          ↓
┌──────────┐  ┌──────────────────┐
│  Loop    │  │ spec_loadfiles() │
│ Ordner & │  │ (einzeln)        │
│ merge    │  └──────────────────┘
└──────────┘           │
      │                │
      └────┬───────────┘
           ↓
┌──────────────────────────┐
│ Alte Threads stoppen     │
│ Nanomap bereinigen       │
│ Frames zerstören         │
│ Garbage Collection       │
└──────────────────────────┘
           │
           ↓
┌──────────────────────────┐
│ os.walk(folder)          │
│ Dateien mit Filter       │
│ sammeln                  │
└──────────────────────────┘
           │
           ↓
┌──────────────────────────┐
│ cmapframe & specframe    │
│ neu erstellen            │
└──────────────────────────┘
           │
           ↓
┌──────────────────────────┐
│ XYMap erstellen mit:     │
│ • files_processed        │
│ • Multi-BG, Linear-BG    │
│ • Cosmic removal options │
│ • Derivate options       │
└──────────────────────────┘
           │
           ↓
┌──────────────────────────┐
│ Power Correction         │
│ (wenn aktiviert)         │
└──────────────────────────┘
           │
           ↓
┌──────────────────────────┐
│ Exporter initialisieren  │
│ → Farbkarte anzeigen     │
│ → Spektren klickbar      │
└──────────────────────────┘
```

### Datei-Sortierung und Merging

```
┌──────────────────────┐
│ specfilesorter       │
│ scan_maindir()       │
└──────────────────────┘
           │
           ↓
┌──────────────────────────────┐
│ os.walk(maindir)             │
│ scan_results =               │
│ [(name, count, path), ...]   │
└──────────────────────────────┘
           │
           ↓
┌──────────────────────────┐
│ Benutzer                 │
│ wählt Ordner aus         │
│ sort_and_process()       │
└──────────────────────────┘
           │
           ↓
┌──────────────────────────┐
│ Für jede Datei:          │
│ shutil.copy2()           │
│ in savedir/name/...      │
│ → _copy_worker()         │
│   Background-Thread      │
└──────────────────────────┘
          │
          ↓
┌──────────────────────────────────┐
│ merge_var == 1  ?                │
│ ├─ Regex-Parse Daten             │
│ ├─ Gruppierung nach prefix+suffix │
│ ├─ Sequenz-Erkennung             │
│ └─ Dateien mergen                │
└──────────────────────────────────┘
           │
           ↓
┌─────────────────┐
│ Fertig!         │
│ Progress -> 100%│
└─────────────────┘
```

---

## Datenformate und Speicherung

### Pickle-Format für XYMap

```python
# Speichern
Nanomap.save_state(filename)
# → filename.pkl (autom. Extension)
# Enthält:
#   ├── Roh-HSI-Daten
#   ├── Verarbeitete Spektren
#   ├── Farbkarte-Bilder
#   ├── ROI-Masken
#   ├── Derivative (falls vorhanden)
#   └── Metadaten

# Laden
Nanomap.load_state(filename)
# → Stellt alle Daten wieder her
# → Ruft update_after_load()
```

### Defaults-Datei

```
defaults.txt
├── data_file (HSI Ordner)
├── filename (Filter)
├── file_extension
├── multiple_Background
├── linear_Background
├── remove_cosmics
├── cosmic_threshold
├── cosmic_width
├── laser_spotsize_nm
├── calculate_first_derivative
├── calculate_second_derivative
├── derivative_polynomial_order
├── derivative_Nfitpoints
├── power_correction
├── loadonstart
└── ... weitere Optionen
```

---

## 🎯 Typische Anwendungsszenarien

### Szenario 1: Einzelne HSI laden und visualisieren
```
1. "Load Data" Tab öffnen
2. Ordner mit Messdaten auswählen
3. Dateiname-Filter eingeben
4. Dateiformat eingeben
5. "Load HSI data" klicken
→ Farbkarte erscheint, Spektren klickbar
```

### Szenario 2: Batch-Verarbeitung über Multiple Tage
```
1. "Process Multiple HSIs" aktivieren
2. Hauptverzeichnis mit Unterlordnern auswählen
3. Speicher-Verzeichnisse setzen
4. "Load HSI data" klicken
→ Für jeden Unterordner:
   - HSI-Bilder generiert
   - Objekte gespeichert (optional)
```

### Szenario 3: Messdaten mit Merging
```
1. "HSI File Sorter" Tab öffnen
2. Hauptverzeichnis auswählen
3. Filter setzen (Dateiname, Extension)
4. "Scan" klicken
5. Resultat übersehen, "Merge consecutive days" aktivieren
6. "Process Selected" klicken
→ Konsekutive Messungen zusammengefügt
```

### Szenario 4: Cosmic Ray & Derivat-Berechnung
```
1. HSI laden (normal)
2. "Remove Cosmics" aktivieren
3. Cosmic-Remove Funktion auswählen
4. Threshold/Width anpassen
5. "calc 1st derivative" und "calc 2nd derivative" aktivieren
6. Polynom-Ordnung wählen
7. "Load HSI data" klicken
→ Alle Verarbeitungen durchgeführt
```

### Szenario 5: Messdaten speichern und später laden
```
# Speichern
1. Nach HSI-Laden
2. "Save the current Data object" nutzen
3. Zielort auswählen → speichern

# Laden
1. "Load a saved Data object"
2. Gespeicherte .pkl-Datei auswählen
→ Kompletter Zustand wiederhergestellt
```

---

## Performance-Überlegungen

### Speicher-Management
- **Garbage Collection:** `gc.collect()` nach Nanomap-Löschung
- **File Handles:** `on_close()` schließt alle offenen Dateien
- **Frame Cleanup:** Alte Frames destruiert vor Neuanlage

### Thread-Sicherheit
- **Copy Worker:** Lock für `_copied_count`
- **GUI Updates:** `after()` für Main-Thread-Sicherheit
- **Stop Events:** Für sauberen Thread-Abbruch

### Large Dataset Handling
- **Canvas Scrolling:** Für 1000+ Spektren
- **Progress Indication:** Benutzer bleibt informiert
- **Hintergrund-Verarbeitung:** Hauptthread bleibt responsiv

---

## Debugging

### Verwendeter Debugger
```python
debugger = DBG.main_Debugger()
```

### Logging
```python
logger.info("Message", extra={'context': 'Context describe'})
```

### Error Engine
```python
error_engine = error_handler.ErrorEngine(
    log_file='logs/specmap.log'
)
```

---

## Zusammenfassung

**main9.py** ist ein ausgefeiltes Tkinter-basiertes Datenverarbeitungs-Tool für hyperspektrale Bilddaten mit folgenden Kernfeatures:

 **Mehrere Datenquellen:** HSI, CLARA, Newton, TCSPC
 **Intelligente Batch-Verarbeitung:** Mehrere HSIs hintereinander
 **Automatisches Merging:** Folgemessungen zusammenfügen
 **Cosmic Ray Entfernung:** Mit verschiedenen Algorithmen
 **Derivat-Berechnung:** Erste und zweite Ableitungen
 **Thread-Management:** Background-Processing ohne UI-Blockierung
 **Persistance:** Speichern/Laden kompletter Verarbeitungszustände
 **Error Handling:** Robustes Error-Handling und Logging
 **Moderne UI:** With Canvas-Scrolling und Mousewheel-Support

---

**Dokumentation erstellt:** 8. April 2026