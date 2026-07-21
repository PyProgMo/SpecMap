import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.colors import ListedColormap
from matplotlib.colors import to_rgba
import tkinter as tk

try:
    from . import deflib1 as deflib
except ImportError:
    import deflib1 as deflib

class Roihandler():
    def __init__(self, roilist={}, pixmatrix=[[]], cmap='viridis', dx=1, dy=1, x0=0, y0=0, update_entries=False, xentry=None, yentry=None):
        self.roi_mode = True
        self.roi_points = []
        self.roi_lines = []
        self.fig = None
        self.roilist = roilist
        self.pixmatrix = pixmatrix
        if cmap not in plt.colormaps():
            cmap = 'viridis' 
        self.cmap = cmap
        self.pixmatrix = np.transpose(self.pixmatrix)
        # selected pixel coordinates
        self.selX = 0
        self.selY = 0
        # x and y axis: 
        self.dx = dx
        self.dy = dy
        self.x0 = x0
        self.y0 = y0
        self.update_entries = update_entries
        self.selectPixX = xentry # tk.Entry widget for selected pixel X coordinate
        self.selectPixY = yentry # tk.Entry widget for selected pixel Y coordinate
        # check if self.selectPixX and self.selectPixY are tk.Entry(xyframe)
        # loading init
        self.construct_roiselgui()
    
    def construct_roiselgui(self, roiselgui=None):
        if roiselgui is not None:
            self.roiselgui = roiselgui

    def construct(self, pixmatrix, roiselgui):
        self.pixmatrix = pixmatrix
        self.pixmatrix = np.transpose(self.pixmatrix)
        self.roiselgui = roiselgui
        self.fig, self.ax = plt.subplots()
        self.fig.subplots_adjust(right=0.89)# distance on right side for buttons
        self.ax.imshow(pixmatrix, cmap=self.cmap)
        # plt.axess([left, bottom, width, height])
        self.ax_button_toggle = plt.axes((0.89, 0.95, 0.1, 0.05))
        self.button_toggle = Button(self.ax_button_toggle, 'Save ROI')
        self.button_toggle.on_clicked(self.toggle_roi)
        self.ax_button_clear = plt.axes((0.89, 0.89, 0.1, 0.05))
        self.button_clear = Button(self.ax_button_clear, 'Clear ROI')
        self.button_clear.on_clicked(self.clear_roi)
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        plt.show()
        self.selnewestroi()

    def toggle_roi(self, event):
        if self.roi_mode == True:
            # Save the current ROI if valid
            if len(self.roi_points) > 2:
                fig, ax = plt.subplots()
                
                # Determine next available ROI index to avoid overwriting if ROIs were deleted
                max_n = 0
                import re
                for k in self.roilist.keys():
                    match = re.match(r'^roi(\d+)', k)
                    if match:
                        n = int(match.group(1))
                        if n > max_n:
                            max_n = n
                for i in range(len(self.roi_points)):
                    self.roi_points[i] = [float(self.roi_points[i][0]), float(self.roi_points[i][1])]
                newroi = deflib.highlight_roi(self.pixmatrix, self.roi_points)
                
                num_pixels = int(np.nansum(newroi))
                new_roi_name = f'roi{max_n + 1} ({num_pixels} px)'
                
                # transpose newroi
                self.roilist[new_roi_name] = newroi
                cax = ax.imshow(newroi, cmap=self.cmap)
                # add colorbar to the plot
                cbar = fig.colorbar(cax, ax=ax)

                plt.show()
                self.roiselgui['values'] = list(self.roilist.keys())
                self.roiselgui.set(new_roi_name)
                # Clear ROI points and lines to allow immediate painting of new ROI
                self.roi_points.clear()
                self.clear_roi_lines()
                # Keep roi_mode = True to allow immediate painting of next ROI
                plt.draw()
            else:
                # Not enough points, just clear and stay in painting mode
                self.roi_points.clear()
                self.clear_roi_lines()
                plt.draw()
        else:
            self.button_toggle.label.set_text('Save ROI')
            self.roi_points.clear()
            self.clear_roi_lines()
            self.roi_mode = True
            plt.draw()

    def clear_roi(self, event):
        self.clear_roi_points()
        self.clear_roi_lines()
        plt.draw()
            
    def on_click(self, event):
        if event.inaxes is None:
            return

        x, y = event.xdata, event.ydata
        if x is None or y is None:
            return

        self.selX = int(round(x))
        self.selY = int(round(y))
        if self.update_entries:
            self.updateselectionentries()

        if self.roi_mode and event.inaxes == self.ax:
            self.roi_points.append((x, y))
            point_plot, = self.ax.plot(x, y, 'ro')
            self.roi_lines.append(point_plot)
            if len(self.roi_points) > 1:
                line_plot, = self.ax.plot([self.roi_points[-2][0], x],
                                            [self.roi_points[-2][1], y], 'r-')
                self.roi_lines.append(line_plot)
            plt.draw()

    def clear_roi_points(self):
        self.roi_points.clear()

    def clear_roi_lines(self):
        for line in self.roi_lines:
            line.remove()
        self.roi_lines.clear()
    
    def plotroi(self, fontsize=12):
        # get selection of self.roiselgui
        roi = self.roilist[self.roiselgui.get()]
        fig, ax = plt.subplots()
        cax = ax.imshow(roi, cmap=self.cmap)
        cbar = fig.colorbar(cax, ax=ax)
        cbar.set_label('ROI', fontsize=fontsize)
        cbar.ax.tick_params(labelsize=fontsize)
        ax.set_title('Region of Interest')
        ax.set_xlabel('Nanostage X Axis in \u03bcm', fontsize=fontsize)
        ax.set_ylabel('Nanostage Y Axis in \u03bcm', fontsize=fontsize)
        if self.update_entries:
            self.bind_click_event(fig, ax)

        plt.show()
    
    def delete_roi(self):
        if self.roiselgui.get() != '':
            del self.roilist[self.roiselgui.get()]
            self.roiselgui['values'] = list(self.roilist.keys())
            self.selnewestroi()
            plt.show()
        else:
            pass
    
    def selnewestroi(self):
        if len(self.roilist) > 0:
            self.roiselgui.set(list(self.roilist.keys())[-1])
        else:
            self.roiselgui.set('')
    
    def on_close(self):
        plt.close('all')

    def _draw_overlay(self, ax, roi, color='red', alpha=0.5):
        """
        Draw ROI as a transparent colored overlay.
        
        Args:
            ax: matplotlib axis
            roi: ROI array
            color: color name (e.g., 'red', 'blue', 'green')
            alpha: transparency level (0-1, default 0.5)
        """
        # Create a masked array to overlay the ROI
        masked_roi = np.ma.masked_where(roi == 0, roi)
        rgba_color = to_rgba(color)
        
        # Create a custom colormap that goes from transparent to the specified color
        # Create colormap: 0 is transparent, any positive value is the color with alpha
        ccmap_colors = [rgba_color]
        ccmap = ListedColormap(ccmap_colors)
        
        # Normalize the ROI values to use the colormap
        ax.imshow(masked_roi, cmap=ccmap, alpha=alpha)

    def _draw_cornerlines(self, ax, roi, color='red'):
        # Outline the exact shape of the ROI connecting its boundaries/corners
        # By providing X and Y meshgrids, we make sure it perfectly aligns with imshow indices
        clean_roi = np.nan_to_num(roi, nan=0)
        Y, X = np.indices(clean_roi.shape)
        ax.contour(X, Y, clean_roi, levels=[0.5], colors=[color], linewidths=3)

    def plot_roi_on_pixmatrix(self, pixmatrix, roiname, vis_type='overlay', color='red', fontsize=12, title=None, unit='A. U.'):
        vis_funcs = {
            'overlay': self._draw_overlay,
            'cornerlines': self._draw_cornerlines
        }

        if roiname in self.roilist:
            roi = self.roilist[roiname]
            fig, ax = plt.subplots()
            
            # Display pixelmatrix as background using standard or passed colormap
            #img = ax.imshow(np.transpose(pixmatrix), cmap='viridis')
            img = ax.imshow(pixmatrix, cmap=self.cmap)
            fig.colorbar(img, ax=ax, label='Intensity')
            
            # Apply the selected visualization function
            if vis_type in vis_funcs:
                vis_funcs[vis_type](ax, roi, color)
            else:
                print(f"Visualization type '{vis_type}' not recognized. Defaulting to 'overlay'.")
                vis_funcs['overlay'](ax, roi, color)
            
            if title is None:
                title = f'Region of Interest: {roiname}'
            ax.set_title(title, fontsize=fontsize)
            ax.set_xlabel('Nanostage X Axis in \u03bcm', fontsize=fontsize)
            ax.set_ylabel('Nanostage Y Axis in \u03bcm', fontsize=fontsize)
            if self.update_entries:
                self.bind_click_event(fig, ax)
            plt.show()
        else:
            print(f"ROI '{roiname}' not found.")
    
    def plot_multiple_rois_on_pixmatrix(self, handler, pixmatrix, roinames, plotmodes, colors, fontsize=14, title=None, unit='A. U.'):
        vis_funcs = {
            'overlay': self._draw_overlay,
            'cornerlines': self._draw_cornerlines
        }
        #display pixmatrix as background using standard colormap before plotting any rois, so it doesn't get overwritten by the rois
        fig, ax = plt.subplots()
        img = ax.imshow(pixmatrix, cmap=self.cmap)
        fig.colorbar(img, ax=ax, label=unit)
        for roiname, plotmode, color in zip(roinames, plotmodes, colors):
            if roiname in handler.roilist:
                roi = handler.roilist[roiname]
                if plotmode in vis_funcs:
                    vis_funcs[plotmode](ax, roi, color)
                else:
                    print(f"Visualization type '{plotmode}' not recognized. Defaulting to 'overlay'.")
                    vis_funcs['overlay'](ax, roi, color)
            else:
                print(f"ROI '{roiname}' not found.")
                
        if title is None:
            title = 'Regions of Interest'
        ax.set_title(title, fontsize=fontsize)
        ax.set_xlabel('Nanostage X Axis in \u03bcm', fontsize=fontsize)
        ax.set_ylabel('Nanostage Y Axis in \u03bcm', fontsize=fontsize)
        if self.update_entries:
            self.bind_click_event(fig, ax) 
        plt.show()

    def bind_click_event(self, fig, ax):
        # click event and update "selected Pixels X and Y"
        fig.canvas.mpl_connect('button_press_event', self.on_click)

    def updateselectionentries(self):
        try:
            self.selectPixX.delete(0, tk.END) #type: ignore
            self.selectPixX.insert(0, str(self.selX)) #type: ignore
            self.selectPixY.delete(0, tk.END) #type: ignore
            self.selectPixY.insert(0, str(self.selY)) #type: ignore
        except Exception as e:
            print(f"Error updating selection entries: {e}. Make sure selectPixX and selectPixY are valid tk.Entry widgets.")

def test_roionpixmatrix():
    pixmatrix = np.random.rand(100, 100)
    
    # Simulate user clicking arbitrary points to form a polygon
    roi_points = [(15, 25), (85, 20), (75, 90), (10, 70)]
    
    # Generate the custom ROI using the real library function
    roi_data = deflib.highlight_roi(np.transpose(pixmatrix), roi_points)
    
    roilist = {'test_roi': roi_data}
    
    handler = Roihandler(roilist=roilist, pixmatrix=pixmatrix)
    
    # Test overlay visualization
    handler.plot_roi_on_pixmatrix(pixmatrix, 'test_roi', vis_type='overlay')
    
    # Test outline (cornerlines) visualization
    handler.plot_roi_on_pixmatrix(pixmatrix, 'test_roi', vis_type='cornerlines', color='red')

def roiindicees2roinames(handler, indicees):
    roinames = []
    # match ROIs by "roiX" substring dynamically
    import re
    for idx in indicees:
        pattern = re.compile(rf"^roi{idx}(?:\s|\(|$)")
        matched = False
        for key in handler.roilist.keys():
            if pattern.match(key):
                roinames.append(key)
                matched = True
                break
        if not matched:
            print(f"Warning: No ROI found matching index '{idx}' (expected 'roi{idx}...')")
    return roinames

if __name__ == "__main__":
    test_roionpixmatrix()