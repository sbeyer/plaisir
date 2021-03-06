/* automatically generated by rust-bindgen 0.59.2 */

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct NodeCoords {
    pub x: f64,
    pub y: f64,
}
#[test]
fn bindgen_test_layout_NodeCoords() {
    assert_eq!(
        ::std::mem::size_of::<NodeCoords>(),
        16usize,
        concat!("Size of: ", stringify!(NodeCoords))
    );
    assert_eq!(
        ::std::mem::align_of::<NodeCoords>(),
        8usize,
        concat!("Alignment of ", stringify!(NodeCoords))
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<NodeCoords>())).x as *const _ as usize },
        0usize,
        concat!(
            "Offset of field: ",
            stringify!(NodeCoords),
            "::",
            stringify!(x)
        )
    );
    assert_eq!(
        unsafe { &(*(::std::ptr::null::<NodeCoords>())).y as *const _ as usize },
        8usize,
        concat!(
            "Offset of field: ",
            stringify!(NodeCoords),
            "::",
            stringify!(y)
        )
    );
}
extern "C" {
    pub fn run(
        dimension: ::std::os::raw::c_int,
        coords: *const NodeCoords,
    ) -> *const ::std::os::raw::c_int;
}
