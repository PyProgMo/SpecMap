import numpy as np
from matplotlib.path import Path

def highlight_roi(Mat, points):
    Mat = np.asarray(Mat)
    result = np.full_like(Mat, np.nan, dtype=float)
    y, x = np.meshgrid(np.arange(Mat.shape[1]), np.arange(Mat.shape[0]))
    points_grid = np.vstack((x.ravel(), y.ravel())).T
    polygon_path = Path(points)
    inside_mask = polygon_path.contains_points(points_grid)
    inside_mask = inside_mask.reshape(Mat.shape)
    result[inside_mask] = 1
    return np.transpose(result)

Mat = np.zeros((100, 100))
points1 = [(10, 10), (10, 20), (20, 20), (20, 10)]
points2 = [(0, 0), (0, 50), (50, 50), (50, 0)]

roi1 = highlight_roi(Mat, points1)
roi2 = highlight_roi(Mat, points2)
print("roi1 sum:", np.nansum(roi1))
print("roi2 sum:", np.nansum(roi2))
