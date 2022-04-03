#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# for colorbar
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

# def par

set_cmp = 'viridis_r'
set_title = 'ricker_wave'

# ============================================================================
def ReadDat(fn, shape=[], dtype="float32" ) -> np.ndarray:
    if ( len(shape) == 0 ):
        return []
        pass
    # if have -1 in shape, read all
    if ( -1 in shape):
        n_elem = -1
    # if no -1 in shape, n_elem = n1*n2*n3....
    else:
        n_elem = np.prod(shape)
        pass  # end if

    arr = np.fromfile(fn, dtype, count=n_elem)
    arr = arr.reshape(shape, order="F")
    return arr

# ============================================================================
def plot_lines():
    nz, nx  = 1024, 512
    data_fn = "/Users/stone/Seismic_code_github/tools/zero_offset_data.dat"
    data    = ReadDat(data_fn, shape=[nz, nx])
    data1   = data[:,255];
    t1 = np.arange(0.0, 2.048, 0.002)
    '''
    t2 = np.arange(0.0, 5.0, 0.02)
    t3 = np.arange(0.0, 2.0, 0.01)
    '''
    box = dict(facecolor='yellow', pad=5, alpha=0.2)
    plt.figure(1)
    # 整个画板的标题
    plt.suptitle('SuperTitle', fontsize=16, fontweight='bold')
    # 位置调整
    plt.subplots_adjust(left=0.2, wspace=0.8, top=0.8)
    # subplot 1
    plt.subplot(1, 1, 1)
    plt.plot(t1, data1, label="$f(x)$")
    plt.title('Subplot 1', color='r')
    plt.ylabel('Y Axis', bbox=box)
    plt.legend()
    '''
    # subplot 2
    plt.subplot(1, 2, 2)
    plt.plot(t3, np.cos(2 * np.pi * t3), '--', label=r"$cos(2\pi x)$")
    plt.title('Subplot 1', color='r')
    plt.xlabel('X Axis', bbox=box)
    plt.ylabel('Y Axis', bbox=box)
    plt.legend()
    '''
    plt.show()

    return


# ============================================================================
def plot_image():

    plt.figure(2)
    ax1 = plt.subplot(1, 1, 1)

    nz, nx  = 251, 81
    data_fn = "/Users/stone/Seismic_code_github/tools/data.dat"
    data    = ReadDat(data_fn, shape=[nz, nx])

    im1 = plt.imshow(data, cmap=set_cmp)
    plt.title(set_title, fontsize=16, fontweight='bold')

    # ------ Add colorbar ------
    # step 1 make ax1 locatable
    ax1_divider = make_axes_locatable(ax1)

    # step 2 split ax1 for a sub ax to create our colorbar
    #       [location]     [width]      [distance to image]
    cax1 = ax1_divider.append_axes("right", size="3%", pad="2%")

    # step 3 call colorbar!
    cb1 = plt.colorbar(im1, cax=cax1)

    plt.show()

    return


# ============================================================================
def plot_multi_sub():
    plt.figure(3)

    plt.subplot(2, 2, 1)
    plt.text(1, 1, "2,2,1")
    plt.xlim([0, 2])
    plt.ylim([0, 2])

    plt.subplot(2, 2, 2)
    plt.text(1, 1, "2,2,2")
    plt.xlim([0, 2])
    plt.ylim([0, 2])

    # using 2nd row as ONE subplot
    plt.subplot(2, 1, 2)
    plt.text(1, 1, "2,1,2")
    plt.xlim([0, 2])
    plt.ylim([0, 2])

    plt.show()


# ============================================================================
def real_main():
    plot_image()
    plot_lines()
    #plot_multi_sub()


# ============================================================================
if __name__ == "__main__":
    real_main()


