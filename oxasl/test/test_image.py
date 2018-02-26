import pytest
import numpy as np

from oxasl.image import AslImage

def test_create_data_singleti():
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage("asldata", data=d, tis=[1.5])
    assert img.ntis == 1
    assert img.tis == [1.5]
    assert not img.plds
    assert img.rpts == [3]
    assert img.order == "prt"

def test_create_data_multiti():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage("asldata", data=d, tis=[1.5, 2.0])
    assert img.ntis == 2
    assert img.tis == [1.5, 2.0]
    assert not img.plds
    assert img.rpts == [2, 2]
    assert img.order == "prt"

def test_create_data_multiti_var_repeats():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage("asldata", data=d, tis=[1.5, 2.0], rpts=[1, 3])
    assert img.ntis == 2
    assert img.tis == [1.5, 2.0]
    assert not img.plds
    assert img.rpts == [1, 3]
    assert img.order == "prt"

def test_create_data_multiphase():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage("asldata", data=d, tis=[1.5], order='mrt', phases=[0, 45, 90, 135, 180, 225, 270, 360])
    assert img.ntis == 1
    assert img.tis == [1.5]
    assert not img.plds
    assert img.rpts == [1]
    assert img.ntc == 8
    assert img.order == "mrt"

def test_create_data_multiphase_nph():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage("asldata", data=d, tis=[1.5], order='mrt', nph=8)
    assert img.ntis == 1
    assert img.tis == [1.5]
    assert not img.plds
    assert img.rpts == [1]
    assert img.ntc == 8
    assert img.order == "mrt"

def test_create_data_diff():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage("asldata", data=d, tis=[1.5], order='rt')
    assert img.ntis == 1
    assert img.tis == [1.5]
    assert not img.plds
    assert img.rpts == [8]
    assert img.ntc == 1
    assert img.order == "rt"

def test_diff_tc():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1], order='prt')
    img = img.diff()
    assert img.ntis == 1
    assert img.tis == [1]
    assert not img.plds
    assert img.rpts == [4]
    assert img.ntc == 1
    assert img.order == "rt"
    data = img.data()
    assert list(data.shape) == [5, 5, 5, 4] 
    assert np.all(data == 1)

def test_diff_ct():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1], order='Prt')
    img = img.diff()
    assert img.ntis == 1
    assert img.tis == [1]
    assert not img.plds
    assert img.rpts == [4]
    assert img.ntc == 1
    assert img.order == "rt"
    data = img.data()
    assert list(data.shape) == [5, 5, 5, 4] 
    assert np.all(data == -1)

def test_reorder_prt_Prt():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1], order='prt')
    img = img.reorder("Prt")
    assert img.ntis == 1
    assert img.tis == [1]
    assert not img.plds
    assert img.rpts == [4]
    assert img.ntc == 2
    assert img.order == "Prt"
    data = img.data()
    assert list(data.shape) == [5, 5, 5, 8] 
    for z in range(8):
        if z % 2 == 0:
            assert np.all(data[:,:,:,z] == z+1)
        else:
            assert np.all(data[:,:,:,z] == z-1)

def test_reorder_prt_rtp():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1], order='prt')
    img = img.reorder("rtp")
    assert img.ntis == 1
    assert img.tis == [1]
    assert not img.plds
    assert img.rpts == [4]
    assert img.ntc == 2
    assert img.order == "rtp"
    data = img.data()
    assert list(data.shape) == [5, 5, 5, 8] 
    for z in range(8):
        if z < 4:
            assert np.all(data[:,:,:,z] == z*2)
        else:
            assert np.all(data[:,:,:,z] == (z-4)*2+1)

def test_reorder_prt_ptr():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1, 2], order='prt')
    img = img.reorder("ptr")
    assert img.ntis == 2
    assert img.tis == [1, 2]
    assert not img.plds
    assert img.rpts == [2, 2]
    assert img.ntc == 2
    assert img.order == "ptr"
    data = img.data()
    assert list(data.shape) == [5, 5, 5, 8] 
    for znew, zold in enumerate([0, 1, 4, 5, 2, 3, 6, 7]):
        assert np.all(data[:,:,:,znew] == zold)

@pytest.mark.xfail
def test_reorder_prt_ptr_var_rpts():
    """ This reordering with variable repeates is not supported yet"""
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1, 2], rpts = [1, 3], order='prt')
    img = img.reorder("ptr")
    assert img.ntis == 2
    assert img.tis == [1, 2]
    assert not img.plds
    assert img.rpts == [1, 3]
    assert img.ntc == 2
    assert img.order == "ptr"
    data = img.data()
    assert list(data.shape) == [5, 5, 5, 8] 
    for znew, zold in enumerate([0, 1, 2, 3, 4, 5, 6, 7]):
        assert np.all(data[:,:,:,znew] == zold)

def test_reorder_prt_rpt():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1, 2], order='prt')
    img = img.reorder("rpt")
    assert img.ntis == 2
    assert img.tis == [1, 2]
    assert not img.plds
    assert img.rpts == [2, 2]
    assert img.ntc == 2
    assert img.order == "rpt"
    data = img.data()
    assert list(data.shape) == [5, 5, 5, 8] 
    for znew, zold in enumerate([0, 2, 1, 3, 4, 6, 5, 7]):
        assert np.all(data[:,:,:,znew] == zold)

def test_reorder_prt_tpr():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1, 2], order='prt')
    img = img.reorder("tpr")
    assert img.ntis == 2
    assert img.tis == [1, 2]
    assert not img.plds
    assert img.rpts == [2, 2]
    assert img.ntc == 2
    assert img.order == "tpr"
    data = img.data()
    assert list(data.shape) == [5, 5, 5, 8] 
    for znew, zold in enumerate([0, 4, 1, 5, 2, 6, 3, 7]):
        assert np.all(data[:,:,:,znew] == zold)

def test_mean_across_repeats_rt():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1, 2], order='rt')
    img = img.mean_across_repeats()
    assert img.ntis == 2
    assert img.tis == [1, 2]
    assert not img.plds
    assert img.rpts == [1, 1]
    assert img.ntc == 1
    assert img.order == "rt"
    data = img.data()
    assert list(data.shape) == [5, 5, 5, 2] 
    for znew, zold in enumerate([1.5, 5.5]):
        assert np.all(data[:,:,:,znew] == zold)

def test_mean_across_repeats_tr():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1, 2], order='tr')
    img = img.mean_across_repeats()
    assert img.ntis == 2
    assert img.tis == [1, 2]
    assert not img.plds
    assert img.rpts == [1, 1]
    assert img.ntc == 1
    assert img.order == "tr"
    data = img.data()
    assert list(data.shape) == [5, 5, 5, 2] 
    for znew, zold in enumerate([3, 4]):
        assert np.all(data[:,:,:,znew] == zold)

def test_split_epochs():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1, 2], order='tr')
    imgs = img.split_epochs(4)
    assert len(imgs) == 2
    for idx, img in enumerate(imgs):
        assert img.ntis == 2
        assert img.tis == [1, 2]
        assert not img.plds
        assert img.rpts == [1, 1]
        assert img.ntc == 1
        assert img.order == "tr"
        data = img.data()
        # Epoch 1 is TIs 1212, data 0123, mean across repeats 12
        # Epoch 2 is TIs 1212, data 4567, mean across repeats 56
        start = idx*4
        for z in range(data.shape[3]):
            assert np.all(data[:,:,:,z] == start+z+1)

def test_split_epochs_overlap():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1, 2], order='tr')
    imgs = img.split_epochs(4, overlap=2)
    assert len(imgs) == 3
    for idx, img in enumerate(imgs):
        assert img.ntis == 2
        assert img.tis == [1, 2]
        assert not img.plds
        assert img.rpts == [1, 1]
        assert img.ntc == 1
        assert img.order == "tr"
        data = img.data()
        # Epoch 1 is TIs 1212, data 0123, mean across repeats 12
        # Epoch 2 is TIs 1212, data 2345, mean across repeats 34
        # Epoch 3 is TIs 1212, data 4567, mean across repeats 56
        start = idx*2
        for z in range(data.shape[3]):
            assert np.all(data[:,:,:,z] == start+z+1)

def test_split_epochs_rt():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1, 2], order='rt')
    imgs = img.split_epochs(4)
    assert len(imgs) == 2
    for idx, img in enumerate(imgs):
        assert img.ntis == 1
        assert img.tis == [1+idx]
        assert not img.plds
        assert img.rpts == [1]
        assert img.ntc == 1
        assert img.order == "rt"
        data = img.data()
        # Epoch 1 is TIs 1111, data 0123, mean across repeats 1.5
        # Epoch 2 is TIs 2222, data 4567, mean across repeats 5.5
        start = idx*4
        for z in range(data.shape[3]):
            assert np.all(data[:,:,:,z] == start+z+1.5)

def test_split_epochs_rt_overlap():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1, 2], order='rt')
    imgs = img.split_epochs(4, overlap=2)
    assert len(imgs) == 3
    for idx, img in enumerate(imgs):
        if idx in (0, 2):
            assert img.ntis == 1
            assert img.tis == [1+idx/2]
            assert not img.plds
            assert img.rpts == [1]
            assert img.ntc == 1
            assert img.order == "rt"
            data = img.data()
            # Epoch 1 is TIs 1111, data 0123, mean across repeats 1.5
            # Epoch 3 is TIs 2222, data 4567, mean across repeats 5.5
            start = idx*2
            for z in range(data.shape[3]):
                assert np.all(data[:,:,:,z] == start+z+1.5)
        else:
            assert img.ntis == 2
            assert img.tis == [1, 2]
            assert not img.plds
            assert img.rpts == [1, 1]
            assert img.ntc == 1
            assert img.order == "rt"
            data = img.data()
            # Epoch 2 is TIs 1122, data 2345, mean across repeats 2.5, 4.5
            for z in range(data.shape[3]):
                assert np.all(data[:,:,:,z] == 2*z+2.5)
            
def test_split_epochs_reorder():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1, 2], order='rt')
    imgs = img.split_epochs(4, time_order="tr")
    assert len(imgs) == 2
    for idx, img in enumerate(imgs):
        assert img.ntis == 2
        assert img.tis == [1, 2]
        assert not img.plds
        assert img.rpts == [1, 1]
        assert img.ntc == 1
        assert img.order == "tr"
        data = img.data()
        # reordered = [0, 4, 1, 5, 2, 6, 3, 7]
        # Epoch 1 is TIs 1212, data 0415, mean across repeats 0.5, 4.5
        # Epoch 2 is TIs 1212, data 2637, mean across repeats 2.5, 6.5
        for z in range(data.shape[3]):
            assert np.all(data[:,:,:,z] == idx*2+0.5+z*4)
            
def test_split_epochs_reorder_overlap():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[:,:,:,z] = z
    img = AslImage("asldata", data=d, tis=[1, 2], order='rt')
    imgs = img.split_epochs(4, overlap=2, time_order="tr")
    assert len(imgs) == 3
    for idx, img in enumerate(imgs):
        assert img.ntis == 2
        assert img.tis == [1, 2]
        assert not img.plds
        assert img.rpts == [1, 1]
        assert img.ntc == 1
        assert img.order == "tr"
        data = img.data()
        # reordered = [0, 4, 1, 5, 2, 6, 3, 7]
        # Epoch 1 is TIs 1212, data 0415, mean across repeats 0.5, 4.5
        # Epoch 2 is TIs 1212, data 1526, mean across repeats 1.5, 5.5
        # Epoch 3 is TIs 1212, data 2637, mean across repeats 2.5, 6.5
        for z in range(data.shape[3]):
            assert np.all(data[:,:,:,z] == idx+0.5+4*z)
