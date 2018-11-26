import pytest
import numpy as np

from fsl.data.image import Image

from oxasl import AslImage

def test_create_data_singleti():
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, tis=[1.5], iaf="tc", order="lrt")
    assert img.ntis == 1
    assert img.tis == [1.5]
    assert not img.have_plds
    assert img.rpts == [3]
    assert img.order == "lrt"

def test_create_data_numtis_single():
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, ntis=1, iaf="tc", order="lrt")
    assert img.ntis == 1
    assert not img.tis
    assert not img.have_plds
    assert img.rpts == [3]
    assert img.order == "lrt"

def test_create_data_multiti():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage(name="asldata", image=d, tis=[1.5, 2.0], iaf="tc", order="lrt")
    assert img.ntis == 2
    assert img.tis == [1.5, 2.0]
    assert not img.have_plds
    assert img.rpts == [2, 2]
    assert img.order == "lrt"

def test_create_data_numtis_multi():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage(name="asldata", image=d, ntis=2, iaf="tc", order="lrt")
    assert img.ntis == 2
    assert not img.tis
    assert not img.have_plds
    assert img.rpts == [2, 2]
    assert img.order == "lrt"

def test_create_data_multiti_var_repeats():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage(name="asldata", image=d, tis=[1.5, 2.0], iaf="tc", order="lrt", rpts=[1, 3])
    assert img.ntis == 2
    assert img.tis == [1.5, 2.0]
    assert not img.have_plds
    assert img.rpts == [1, 3]
    assert img.order == "lrt"

def test_create_data_multiphase():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage(name="asldata", image=d, tis=[1.5], iaf='mp', order='lrt', phases=[0, 45, 90, 135, 180, 225, 270, 315])
    assert img.ntis == 1
    assert img.tis == [1.5]
    assert not img.have_plds
    assert img.rpts == [1]
    assert img.ntc == 8
    assert img.phases == [0, 45, 90, 135, 180, 225, 270, 315]
    assert img.order == "lrt"

def test_create_data_multiphase_nphases():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage(name="asldata", image=d, tis=[1.5], iaf="mp", order='lrt', nphases=8)
    assert img.ntis == 1
    assert img.tis == [1.5]
    assert not img.have_plds
    assert img.rpts == [1]
    assert img.ntc == 8
    assert img.phases == [0, 45, 90, 135, 180, 225, 270, 315]
    assert img.order == "lrt"

def test_create_data_plds():
    d = np.random.rand(5, 5, 5, 6)
    img = AslImage(name="asldata", image=d, plds=[1.5], tau=1.8, iaf="tc", order="lrt")
    assert img.ntis == 1
    assert img.plds == [1.5]
    assert img.tis == [pytest.approx(3.3)]
    assert img.have_plds
    assert img.casl
    assert img.rpts == [3]
    assert img.order == "lrt"

def test_create_data_multiplds():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage(name="asldata", image=d, plds=[1.5, 1.6], tau=1.8, iaf="tc", order="lrt")
    assert img.ntis == 2
    assert img.plds == [1.5, 1.6]
    assert img.tis == [pytest.approx(3.3), pytest.approx(3.4)]
    assert img.have_plds
    assert img.casl
    assert img.rpts == [2, 2]
    assert img.order == "lrt"

def test_create_data_multitaus():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage(name="asldata", image=d, plds=[1.5, 1.6], taus=[1.8, 2.0], iaf="tc", order="lrt")
    assert img.ntis == 2
    assert img.plds == [1.5, 1.6]
    assert img.tis == [pytest.approx(3.3), pytest.approx(3.6)]
    assert img.have_plds
    assert img.casl
    assert img.rpts == [2, 2]
    assert img.order == "lrt"

def test_create_data_strtaus():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage(name="asldata", image=d, plds=[1.5, 1.6], taus="1.8, 2.0", iaf="tc", order="lrt")
    assert img.ntis == 2
    assert img.plds == [1.5, 1.6]
    assert img.tis == [pytest.approx(3.3), pytest.approx(3.6)]
    assert img.have_plds
    assert img.casl
    assert img.rpts == [2, 2]
    assert img.order == "lrt"

def test_create_data_floattaus():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage(name="asldata", image=d, plds=[1.5, 1.6], tau=1.8, iaf="tc", order="lrt")
    assert img.ntis == 2
    assert img.plds == [1.5, 1.6]
    assert img.tis == [pytest.approx(3.3), pytest.approx(3.4)]
    assert img.have_plds
    assert img.casl
    assert img.rpts == [2, 2]
    assert img.order == "lrt"

def test_create_data_inttaus():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage(name="asldata", image=d, plds=[1.5, 1.6], tau=1, iaf="tc", order="lrt")
    assert img.ntis == 2
    assert img.plds == [1.5, 1.6]
    assert img.tis == [pytest.approx(2.5), pytest.approx(2.6)]
    assert img.have_plds
    assert img.casl
    assert img.rpts == [2, 2]
    assert img.order == "lrt"

def test_create_data_plds_pasl():
    """ Odd usage, PLDs treated as TIs in non-CASL acquisitions """
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage(name="asldata", image=d, plds=[1.5, 2.5], tau=1.8, casl=False, iaf="tc", order="lrt")
    assert img.ntis == 2
    assert img.tis == [1.5, 2.5]
    assert img.plds == [1.5, 2.5]
    assert img.have_plds
    assert img.rpts == [2, 2]
    assert img.order == "lrt"

def test_create_data_diff():
    d = np.random.rand(5, 5, 5, 8)
    img = AslImage(name="asldata", image=d, tis=[1.5], iaf="diff", order='rt')
    assert img.ntis == 1
    assert img.tis == [1.5]
    assert not img.have_plds
    assert img.rpts == [8]
    assert img.ntc == 1
    assert img.order == "rt"

def test_diff_tc():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1], iaf="tc", order='lrt')
    img = img.diff()
    assert img.ntis == 1
    assert img.tis == [1]
    assert not img.have_plds
    assert img.rpts == [4]
    assert img.ntc == 1
    assert img.order == "rt"
    data = img.nibImage.get_data()
    assert list(data.shape) == [5, 5, 5, 4] 
    assert np.all(data == 1)

def test_diff_ct():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1], iaf="ct", order='lrt')
    img = img.diff()
    assert img.ntis == 1
    assert img.tis == [1]
    assert not img.have_plds
    assert img.rpts == [4]
    assert img.ntc == 1
    assert img.order == "rt"
    data = img.nibImage.get_data()
    assert list(data.shape) == [5, 5, 5, 4] 
    assert np.all(data == -1)

def test_reorder_tc_ct():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1], iaf="tc", order='lrt')
    img = img.reorder(iaf="ct")
    assert img.ntis == 1
    assert img.tis == [1]
    assert not img.have_plds
    assert img.rpts == [4]
    assert img.ntc == 2
    assert img.order == "lrt"
    assert img.iaf == "ct"
    data = img.nibImage.get_data()
    assert list(data.shape) == [5, 5, 5, 8] 
    for z in range(8):
        if z % 2 == 0:
            assert np.all(data[..., z] == z+1)
        else:
            assert np.all(data[..., z] == z-1)

def test_reorder_lrt_rtl():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1], iaf="tc", order='lrt')
    img = img.reorder("rtl")
    assert img.ntis == 1
    assert img.tis == [1]
    assert not img.have_plds
    assert img.rpts == [4]
    assert img.ntc == 2
    assert img.order == "rtl"
    data = img.nibImage.get_data()
    assert list(data.shape) == [5, 5, 5, 8] 
    for z in range(8):
        if z < 4:
            assert np.all(data[..., z] == z*2)
        else:
            assert np.all(data[..., z] == (z-4)*2+1)

def test_reorder_lrt_ltr():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1, 2], iaf="tc", order='lrt')
    img = img.reorder("ltr")
    assert img.ntis == 2
    assert img.tis == [1, 2]
    assert not img.have_plds
    assert img.rpts == [2, 2]
    assert img.ntc == 2
    assert img.order == "ltr"
    data = img.nibImage.get_data()
    assert list(data.shape) == [5, 5, 5, 8] 
    for znew, zold in enumerate([0, 1, 4, 5, 2, 3, 6, 7]):
        assert np.all(data[..., znew] == zold)

def test_reorder_lrt_ltr_var_rpts():
    """ 
    This reordering with variable repeats is not supported. In priciple
    it is possible (TI1_R1, TI2_R1, TI2_R2, TI2_R3) but this seems unlikely
    and is probably more likely an error
    """
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1, 2], rpts = [1, 3], iaf="tc", order='lrt')
    with pytest.raises(Exception):
        img.reorder("ltr")
    #assert img.ntis == 2
    #assert img.tis == [1, 2]
    #assert not img.have_plds
    #assert img.rpts == [1, 3]
    #assert img.ntc == 2
    #assert img.order == "ltr"
    #data = img.nibImage.get_data()
    #assert list(data.shape) == [5, 5, 5, 8] 
    #for znew, zold in enumerate([0, 1, 2, 3, 4, 5, 6, 7]):
    #    assert np.all(data[..., znew] == zold)

def test_reorder_lrt_rlt():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1, 2], iaf="tc", order='lrt')
    img = img.reorder("rlt")
    assert img.ntis == 2
    assert img.tis == [1, 2]
    assert not img.have_plds
    assert img.rpts == [2, 2]
    assert img.ntc == 2
    assert img.order == "rlt"
    data = img.nibImage.get_data()
    assert list(data.shape) == [5, 5, 5, 8] 
    for znew, zold in enumerate([0, 2, 1, 3, 4, 6, 5, 7]):
        assert np.all(data[..., znew] == zold)

def test_reorder_lrt_tlr():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1, 2], iaf="tc", order='lrt')
    img = img.reorder("tlr")
    assert img.ntis == 2
    assert img.tis == [1, 2]
    assert not img.have_plds
    assert img.rpts == [2, 2]
    assert img.ntc == 2
    assert img.order == "tlr"
    data = img.nibImage.get_data()
    assert list(data.shape) == [5, 5, 5, 8] 
    for znew, zold in enumerate([0, 4, 1, 5, 2, 6, 3, 7]):
        assert np.all(data[..., znew] == zold)

def test_mean_across_repeats_rt():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1, 2], order='rt')
    img = img.mean_across_repeats()
    assert img.ntis == 2
    assert img.tis == [1, 2]
    assert not img.have_plds
    assert img.rpts == [1, 1]
    assert img.ntc == 1
    assert img.order == "rt"
    data = img.nibImage.get_data()
    assert list(data.shape) == [5, 5, 5, 2] 
    for znew, zold in enumerate([1.5, 5.5]):
        assert np.all(data[..., znew] == zold)

def test_mean_across_repeats_tr():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1, 2], order='tr')
    img = img.mean_across_repeats()
    assert img.ntis == 2
    assert img.tis == [1, 2]
    assert not img.have_plds
    assert img.rpts == [1, 1]
    assert img.ntc == 1
    assert img.order == "tr"
    data = img.nibImage.get_data()
    assert list(data.shape) == [5, 5, 5, 2] 
    for znew, zold in enumerate([3, 4]):
        assert np.all(data[..., znew] == zold)

def test_mean_across_repeats_var_rpts():
    d = np.random.rand(5, 5, 5, 86)
    #for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1, 2, 3, 4, 5], rpts=[6, 6, 6, 10, 15], iaf="tc", order='lrt')
    img = img.mean_across_repeats()
    assert img.ntis == 5
    assert img.tis == [1, 2, 3, 4, 5]
    assert img.rpts == [1, 1, 1, 1, 1]
    assert img.ntc == 1
    assert img.order == "rt"
    data = img.nibImage.get_data()
    assert list(data.shape) == [5, 5, 5, 5] 
    #for znew, zold in enumerate([3, 4]):
    #    assert np.all(data[..., znew] == zold)

def test_perf_weighted_tr():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1, 2], order='tr')
    pwi = img.perf_weighted()
    assert(isinstance(pwi, Image))
    data = pwi.data
    assert list(data.shape) == [5, 5, 5] 
    assert np.all(np.mean(d, -1) == data)

def test_vol_idx_tr():
    d = np.zeros([5, 5, 5, 8])
    img = AslImage(name="asldata", image=d, tis=[1, 2], order='tr')
    assert img.get_vol_index(0, 0, 0) == 0
    assert img.get_vol_index(0, 1, 0) == 1
    assert img.get_vol_index(0, 0, 1) == 2
    assert img.get_vol_index(0, 1, 1) == 3
    assert img.get_vol_index(0, 0, 2) == 4
    assert img.get_vol_index(0, 1, 2) == 5
    assert img.get_vol_index(0, 0, 3) == 6
    assert img.get_vol_index(0, 1, 3) == 7
    
def test_vol_idx_rt():
    d = np.zeros([5, 5, 5, 8])
    img = AslImage(name="asldata", image=d, tis=[1, 2], order='rt')
    assert img.get_vol_index(0, 0, 0) == 0
    assert img.get_vol_index(0, 1, 0) == 4
    assert img.get_vol_index(0, 0, 1) == 1
    assert img.get_vol_index(0, 1, 1) == 5
    assert img.get_vol_index(0, 0, 2) == 2
    assert img.get_vol_index(0, 1, 2) == 6
    assert img.get_vol_index(0, 0, 3) == 3
    assert img.get_vol_index(0, 1, 3) == 7
    
def test_vol_idx_lrt():
    d = np.zeros([5, 5, 5, 8])
    img = AslImage(name="asldata", image=d, tis=[1, 2], iaf="tc", order='lrt')
    assert img.get_vol_index(0, 0, 0) == 0
    assert img.get_vol_index(1, 0, 0) == 1
    assert img.get_vol_index(0, 1, 0) == 4
    assert img.get_vol_index(1, 1, 0) == 5
    assert img.get_vol_index(0, 0, 1) == 2
    assert img.get_vol_index(1, 0, 1) == 3
    assert img.get_vol_index(0, 1, 1) == 6
    assert img.get_vol_index(1, 1, 1) == 7
    
def test_vol_idx_ltr():
    d = np.zeros([5, 5, 5, 8])
    img = AslImage(name="asldata", image=d, tis=[1, 2], iaf="tc", order='ltr')
    assert img.get_vol_index(0, 0, 0) == 0
    assert img.get_vol_index(1, 0, 0) == 1
    assert img.get_vol_index(0, 1, 0) == 2
    assert img.get_vol_index(1, 1, 0) == 3
    assert img.get_vol_index(0, 0, 1) == 4
    assert img.get_vol_index(1, 0, 1) == 5
    assert img.get_vol_index(0, 1, 1) == 6
    assert img.get_vol_index(1, 1, 1) == 7
       
def test_vol_idx_var_rpts_rt():
    d = np.zeros([5, 5, 5, 7])
    img = AslImage(name="asldata", image=d, tis=[1, 2], rpts=[3, 4], iaf="diff", order='rt')
    assert img.get_vol_index(0, 0, 0) == 0
    assert img.get_vol_index(0, 0, 1) == 1
    assert img.get_vol_index(0, 0, 2) == 2
    assert img.get_vol_index(0, 1, 0) == 3
    assert img.get_vol_index(0, 1, 1) == 4
    assert img.get_vol_index(0, 1, 2) == 5
    assert img.get_vol_index(0, 1, 3) == 6
       
def test_vol_idx_var_rpts_tr():
    d = np.zeros([5, 5, 5, 7])
    img = AslImage(name="asldata", image=d, tis=[1, 2], rpts=[3, 4], iaf="diff", order='tr')
    assert img.get_vol_index(0, 0, 0) == 0
    assert img.get_vol_index(0, 1, 0) == 1
    assert img.get_vol_index(0, 0, 1) == 2
    assert img.get_vol_index(0, 1, 1) == 3
    assert img.get_vol_index(0, 0, 2) == 4
    assert img.get_vol_index(0, 1, 2) == 5
    assert img.get_vol_index(0, 1, 3) == 6
       
def test_vol_idx_var_rpts_lrt():
    d = np.zeros([5, 5, 5, 10])
    img = AslImage(name="asldata", image=d, tis=[1, 2], rpts=[3, 2], iaf="tc", order='lrt')
    assert img.get_vol_index(0, 0, 0) == 0
    assert img.get_vol_index(1, 0, 0) == 1
    assert img.get_vol_index(0, 0, 1) == 2
    assert img.get_vol_index(1, 0, 1) == 3
    assert img.get_vol_index(0, 0, 2) == 4
    assert img.get_vol_index(1, 0, 2) == 5
    assert img.get_vol_index(0, 1, 0) == 6
    assert img.get_vol_index(1, 1, 0) == 7
    assert img.get_vol_index(0, 1, 1) == 8
    assert img.get_vol_index(1, 1, 1) == 9
    
def test_split_epochs():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1, 2], order='tr')
    imgs = img.split_epochs(4)
    assert len(imgs) == 2
    for idx, img in enumerate(imgs):
        assert img.ntis == 2
        assert img.tis == [1, 2]
        assert not img.have_plds
        assert img.rpts == [1, 1]
        assert img.ntc == 1
        assert img.order == "tr"
        data = img.nibImage.get_data()
        # Epoch 1 is TIs 1212, data 0123, mean across repeats 12
        # Epoch 2 is TIs 1212, data 4567, mean across repeats 56
        start = idx*4
        for z in range(data.shape[3]):
            assert np.all(data[..., z] == start+z+1)

def test_split_epochs_overlap():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1, 2], order='tr')
    imgs = img.split_epochs(4, overlap=2)
    assert len(imgs) == 3
    for idx, img in enumerate(imgs):
        assert img.ntis == 2
        assert img.tis == [1, 2]
        assert not img.have_plds
        assert img.rpts == [1, 1]
        assert img.ntc == 1
        assert img.order == "tr"
        data = img.nibImage.get_data()
        # Epoch 1 is TIs 1212, data 0123, mean across repeats 12
        # Epoch 2 is TIs 1212, data 2345, mean across repeats 34
        # Epoch 3 is TIs 1212, data 4567, mean across repeats 56
        start = idx*2
        for z in range(data.shape[3]):
            assert np.all(data[..., z] == start+z+1)

def test_split_epochs_rt():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1, 2], order='rt')
    imgs = img.split_epochs(4)
    assert len(imgs) == 2
    for idx, img in enumerate(imgs):
        assert img.ntis == 1
        assert img.tis == [1+idx]
        assert not img.have_plds
        assert img.rpts == [1]
        assert img.ntc == 1
        assert img.order == "rt"
        data = img.nibImage.get_data()
        # Epoch 1 is TIs 1111, data 0123, mean across repeats 1.5
        # Epoch 2 is TIs 2222, data 4567, mean across repeats 5.5
        start = idx*4
        for z in range(data.shape[3]):
            assert np.all(data[..., z] == start+z+1.5)

def test_split_epochs_rt_overlap():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1, 2], order='rt')
    imgs = img.split_epochs(4, overlap=2)
    assert len(imgs) == 3
    for idx, img in enumerate(imgs):
        if idx in (0, 2):
            assert img.ntis == 1
            assert img.tis == [1+idx/2]
            assert not img.have_plds
            assert img.rpts == [1]
            assert img.ntc == 1
            assert img.order == "rt"
            data = img.nibImage.get_data()
            # Epoch 1 is TIs 1111, data 0123, mean across repeats 1.5
            # Epoch 3 is TIs 2222, data 4567, mean across repeats 5.5
            start = idx*2
            for z in range(data.shape[3]):
                assert np.all(data[..., z] == start+z+1.5)
        else:
            assert img.ntis == 2
            assert img.tis == [1, 2]
            assert not img.have_plds
            assert img.rpts == [1, 1]
            assert img.ntc == 1
            assert img.order == "rt"
            data = img.nibImage.get_data()
            # Epoch 2 is TIs 1122, data 2345, mean across repeats 2.5, 4.5
            for z in range(data.shape[3]):
                assert np.all(data[..., z] == 2*z+2.5)
            
def test_split_epochs_reorder():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1, 2], order='rt')
    imgs = img.split_epochs(4, time_order="tr")
    assert len(imgs) == 2
    for idx, img in enumerate(imgs):
        assert img.ntis == 2
        assert img.tis == [1, 2]
        assert not img.have_plds
        assert img.rpts == [1, 1]
        assert img.ntc == 1
        assert img.order == "tr"
        data = img.nibImage.get_data()
        # reordered = [0, 4, 1, 5, 2, 6, 3, 7]
        # Epoch 1 is TIs 1212, data 0415, mean across repeats 0.5, 4.5
        # Epoch 2 is TIs 1212, data 2637, mean across repeats 2.5, 6.5
        for z in range(data.shape[3]):
            assert np.all(data[..., z] == idx*2+0.5+z*4)
            
def test_split_epochs_reorder_overlap():
    d = np.zeros([5, 5, 5, 8])
    for z in range(8): d[..., z] = z
    img = AslImage(name="asldata", image=d, tis=[1, 2], order='rt')
    imgs = img.split_epochs(4, overlap=2, time_order="tr")
    assert len(imgs) == 3
    for idx, img in enumerate(imgs):
        assert img.ntis == 2
        assert img.tis == [1, 2]
        assert not img.have_plds
        assert img.rpts == [1, 1]
        assert img.ntc == 1
        assert img.order == "tr"
        data = img.nibImage.get_data()
        # reordered = [0, 4, 1, 5, 2, 6, 3, 7]
        # Epoch 1 is TIs 1212, data 0415, mean across repeats 0.5, 4.5
        # Epoch 2 is TIs 1212, data 1526, mean across repeats 1.5, 5.5
        # Epoch 3 is TIs 1212, data 2637, mean across repeats 2.5, 6.5
        for z in range(data.shape[3]):
            assert np.all(data[..., z] == idx+0.5+4*z)

def test_derived():
    d1 = np.random.rand(5, 5, 5, 8)
    d2 = np.random.rand(5, 5, 5, 8)
    img = AslImage(name="d1", image=d1, plds=[1, 2], iaf="ct", ibf="tis", casl=True, bolus=[3, 4])
    assert img.ntis == 2
    assert img.have_plds
    assert img.plds == [1, 2]
    assert img.tis == [4, 6]
    assert img.rpts == [2, 2]
    assert img.ntc == 2
    assert img.order == "lrt"
    
    img2 = img.derived(d2, name="d2")
    assert img2.ntis == 2
    assert img2.have_plds
    assert img2.plds == [1, 2]
    assert img2.tis == [4, 6]
    assert img2.rpts == [2, 2]
    assert img2.ntc == 2
    assert img2.order == "lrt"
    