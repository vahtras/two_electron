use numpy::{PyArray1, PyArray2, PyArray3, PyArrayMethods, PyUntypedArrayMethods};
use pyo3::prelude::*;

#[pyfunction]
fn fck(
    j: &Bound<'_, PyArray2<f64>>,
    k: &Bound<'_, PyArray2<f64>>,
    dj: &Bound<'_, PyArray2<f64>>,
    dk: &Bound<'_, PyArray2<f64>>,
    buf: &Bound<'_, PyArray1<f64>>,
    ibuf: &Bound<'_, PyArray2<i64>>
    ) {


    let buf = buf.readonly();
    let buf = buf.as_array();

    let ibuf = ibuf.readonly();
    let ibuf = ibuf.as_array();

    let dj = dj.readonly();
    let dj = dj.as_array();

    let dk = dk.readonly();
    let dk = dk.as_array();

    let mut j = j.readwrite();
    let mut j = j.as_array_mut();
    let mut k = k.readwrite();
    let mut k = k.as_array_mut();

    for (index, gint) in buf.iter().enumerate() {
        let s = (ibuf[[0, index]] - 1) as usize;
        let r = (ibuf[[1, index]] - 1) as usize;
        let q = (ibuf[[2, index]] - 1) as usize;
        let p = (ibuf[[3, index]] - 1) as usize;
        let mut g = *gint;

        if p == q { g *= 0.5};
        if r == s { g *= 0.5};
        if p == r && s == q { g *= 0.5 };

        let pq: [usize; 2] = [p, q];
        let qp: [usize; 2] = [q, p];
        let rs: [usize; 2] = [r, s];
        let sr: [usize; 2] = [s, r];


        let fadd = g*(dj[rs] + dj[rs]);
        j[pq] += fadd;
        j[qp] += fadd;
        let fadd = g*(dj[pq] + dj[pq]);
        j[rs] += fadd;
        j[sr] += fadd;

        let ps: [usize; 2] = [p, s];
        let sp: [usize; 2] = [s, p];
        let pr: [usize; 2] = [p, r];
        let rp: [usize; 2] = [r, p];
        let qr: [usize; 2] = [q, r];
        let rq: [usize; 2] = [r, q];
        let qs: [usize; 2] = [q, s];
        let sq: [usize; 2] = [s, q];

        k[ps] += g*dk[rq];
        k[pr] += g*dk[sq];
        k[qs] += g*dk[rp];
        k[qr] += g*dk[sp];
        k[rq] += g*dk[ps];
        k[sq] += g*dk[pr];
        k[rp] += g*dk[qs];
        k[sp] += g*dk[qr];
    }


}

#[pyfunction]
fn fckab(
    j: &Bound<'_, PyArray3<f64>>,
    ka: &Bound<'_, PyArray3<f64>>,
    kb: &Bound<'_, PyArray3<f64>>,
    da: &Bound<'_, PyArray3<f64>>,
    db: &Bound<'_, PyArray3<f64>>,
    buf: &Bound<'_, PyArray1<f64>>,
    ibuf: &Bound<'_, PyArray2<i64>>
    ) {

    println!("{:?}", buf.len());
    println!("{:?}", ibuf.shape());
    println!("{:?}", j.shape());
    let nd = j.shape()[2];

    let buf = buf.readonly();
    let buf = buf.as_array();

    let ibuf = ibuf.readonly();
    let ibuf = ibuf.as_array();

    let da = da.readonly();
    let da = da.as_array();

    let db = db.readonly();
    let db = db.as_array();

    let mut j = j.readwrite();
    let mut j = j.as_array_mut();
    let mut ka = ka.readwrite();
    let mut ka = ka.as_array_mut();
    let mut kb = kb.readwrite();
    let mut kb = kb.as_array_mut();

    for (index, gint) in buf.iter().enumerate() {
        let s = (ibuf[[0, index]] - 1) as usize;
        let r = (ibuf[[1, index]] - 1) as usize;
        let q = (ibuf[[2, index]] - 1) as usize;
        let p = (ibuf[[3, index]] - 1) as usize;
        let mut g = *gint;
        println!("({p}{q}|{r}{s})={g}");
        if p == q { g *= 0.5};
        if r == s { g *= 0.5};
        if p == r && s == q { g *= 0.5 };

        for i in 0..nd {
            let pq: [usize; 3] = [p, q, i];
            let qp: [usize; 3] = [q, p, i];
            let rs: [usize; 3] = [r, s, i];
            let sr: [usize; 3] = [s, r, i];


            let fadd = g*(da[rs] + db[rs] + da[sr] + db[sr]);
            j[pq] += fadd;
            j[qp] += fadd;
            let fadd = g*(da[pq] + db[pq] + da[qp] + db[qp]);
            j[rs] += fadd;
            j[sr] += fadd;

            let ps: [usize; 3] = [p, s, i];
            let sp: [usize; 3] = [s, p, i];
            let pr: [usize; 3] = [p, r, i];
            let rp: [usize; 3] = [r, p, i];
            let qr: [usize; 3] = [q, r, i];
            let rq: [usize; 3] = [r, q, i];
            let qs: [usize; 3] = [q, s, i];
            let sq: [usize; 3] = [s, q, i];

            ka[ps] += g*da[rq];
            ka[pr] += g*da[sq];
            ka[qs] += g*da[rp];
            ka[qr] += g*da[sp];
            ka[rq] += g*da[ps];
            ka[sq] += g*da[pr];
            ka[rp] += g*da[qs];
            ka[sp] += g*da[qr];

            kb[ps] += g*db[rq];
            kb[pr] += g*db[sq];
            kb[qs] += g*db[rp];
            kb[qr] += g*db[sp];
            kb[rq] += g*db[ps];
            kb[sq] += g*db[pr];
            kb[rp] += g*db[qs];
            kb[sp] += g*db[qr];

        }
    }


}

#[pymodule]
fn sirfck(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(fck, m)?)?;
    m.add_function(wrap_pyfunction!(fckab, m)?)?;
    Ok(())
}
