// fn calculate_weight_area(
//     &self,
//     ptarget: &Point,
//     edges: (usize, usize, usize),
// ) -> Result<f64, InterpolatorError> {
//     let point_prev = &self.points[self.triangles[edges.0]];
//     let point_base = &self.points[self.triangles[edges.1]];
//     let point_next = &self.points[self.triangles[edges.2]];

//     let mprev = &Point {
//         x: (point_base.x + point_prev.x) / 2.,
//         y: (point_base.y + point_prev.y) / 2.,
//     };
//     let mnext = &Point {
//         x: (point_base.x + point_next.x) / 2.,
//         y: (point_base.y + point_next.y) / 2.,
//     };

//     let mut ce = edges.0;

//     let pre = {
//         let mut pre = 0.;
//         let mut cs1 = mprev.clone();
//         for dcount in 0..self.degree_limitation {
//             let cit = ce / 3;
//             let triangle = [
//                 &self.points[self.triangles[cit * 3]],
//                 &self.points[self.triangles[cit * 3 + 1]],
//                 &self.points[self.triangles[cit * 3 + 2]],
//             ];
//             let c = circumcenter(&triangle);
//             pre += (cs1.x - c.x) * (cs1.y + c.y);
//             cs1 = c;
//             let next = next_harfedge(ce);
//             if edges.1 == next {
//                 break;
//             }
//             ce = self.harfedges[next];

//             if self.detect_too_large_degree(dcount) {
//                 return Err(InterpolatorError::TooManyNeighbors(self.degree_limitation));
//             }
//         }
//         pre + (cs1.x - mnext.x) * (cs1.y + mnext.y) + (mnext.x - mprev.x) * (mnext.y + mprev.y)
//     };

//     let gprev = circumcenter(&[ptarget, point_base, point_prev]);
//     let gnext = circumcenter(&[ptarget, point_base, point_next]);

//     let post = (mprev.x - gprev.x) * (mprev.y + gprev.y)
//         + (gprev.x - gnext.x) * (gprev.y + gnext.y)
//         + (gnext.x - mnext.x) * (gnext.y + mnext.y)
//         + (mnext.x - mprev.x) * (mnext.y + mprev.y);

//     Ok(pre - post)
// }
