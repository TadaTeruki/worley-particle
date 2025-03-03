use std::error::Error;

use contour::ContourBuilder;

use crate::map::ParticleMap;

use super::{InterpolationMethod, ParticleMapAttributeLerp};

pub struct Band {
    pub threshold: f64,
    pub polygons: Vec<Vec<(f64, f64)>>,
}

enum VectorizationType {
    Contour,
    Isoband,
}

impl<T: ParticleMapAttributeLerp + Into<f64>> ParticleMap<T> {
    fn vectorize(
        &self,
        vectorization_type: VectorizationType,
        corners: ((f64, f64), (f64, f64)),
        rasterise_scale: f64,
        thresholds: &[f64],
        interp_method: &InterpolationMethod,
        smooth: bool,
    ) -> Result<Vec<Band>, Box<dyn Error>> {
        let ((original_min_x, original_min_y), (original_max_x, original_max_y)) = corners;

        let expansion = rasterise_scale * self.params.scale;

        let raster_min_x = (original_min_x * expansion).floor() as i64;
        let raster_max_x = (original_max_x * expansion).ceil() as i64;
        let raster_min_y = (original_min_y * expansion).floor() as i64;
        let raster_max_y = (original_max_y * expansion).ceil() as i64;

        let domain_min_x = raster_min_x as f64 / expansion;
        let domain_max_x = raster_max_x as f64 / expansion;
        let domain_min_y = raster_min_y as f64 / expansion;
        let domain_max_y = raster_max_y as f64 / expansion;

        let width = (raster_max_x - raster_min_x) as usize;
        let height = (raster_max_y - raster_min_y) as usize;
        let raster = self.rasterise(
            width,
            height,
            ((domain_min_x, domain_min_y), (domain_max_x, domain_max_y)),
            interp_method,
        );

        let raster_f64 = raster
            .iter()
            .flat_map(|row| {
                row.iter()
                    .map(|value: &Option<T>| {
                        if let Some(value) = value {
                            value.clone().into()
                        } else {
                            f64::MIN
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        let result_coord_to_domain = |x: f64, y: f64| {
            (
                (x / width as f64) * (domain_max_x - domain_min_x) + domain_min_x,
                (y / height as f64) * (domain_max_y - domain_min_y) + domain_min_y,
            )
        };

        match vectorization_type {
            VectorizationType::Contour => {
                let result =
                    ContourBuilder::new(width, height, smooth).contours(&raster_f64, thresholds);
                if let Err(e) = result {
                    return Err(e.to_string().into());
                }

                let bands = result
                    .unwrap()
                    .iter()
                    .map(|contour| {
                        let polygons = contour
                            .geometry()
                            .iter()
                            .map(|polygon| {
                                let exterior = polygon.exterior();
                                exterior
                                    .0
                                    .iter()
                                    .map(|coord| result_coord_to_domain(coord.x, coord.y))
                                    .collect::<Vec<_>>()
                            })
                            .collect::<Vec<_>>();

                        Band {
                            threshold: contour.threshold(),
                            polygons,
                        }
                    })
                    .collect::<Vec<_>>();

                Ok(bands)
            }
            VectorizationType::Isoband => {
                let result =
                    ContourBuilder::new(width, height, smooth).isobands(&raster_f64, thresholds);
                if let Err(e) = result {
                    return Err(e.to_string().into());
                }
                let bands = result
                    .unwrap()
                    .iter()
                    .map(|band| {
                        let polygons = band
                            .geometry()
                            .iter()
                            .map(|polygon| {
                                let exterior = polygon.exterior();
                                exterior
                                    .0
                                    .iter()
                                    .map(|coord| result_coord_to_domain(coord.x, coord.y))
                                    .collect::<Vec<_>>()
                            })
                            .collect::<Vec<_>>();

                        Band {
                            threshold: band.min_v(),
                            polygons,
                        }
                    })
                    .collect::<Vec<_>>();

                Ok(bands)
            }
        }
    }

    pub fn contours(
        &self,
        corners: ((f64, f64), (f64, f64)),
        rasterise_scale: f64,
        thresholds: &[f64],
        interp_method: &InterpolationMethod,
        smooth: bool,
    ) -> Result<Vec<Band>, Box<dyn Error>> {
        self.vectorize(
            VectorizationType::Contour,
            corners,
            rasterise_scale,
            thresholds,
            interp_method,
            smooth,
        )
    }

    pub fn isobands(
        &self,
        corners: ((f64, f64), (f64, f64)),
        rasterise_scale: f64,
        thresholds: &[f64],
        interp_method: &InterpolationMethod,
        smooth: bool,
    ) -> Result<Vec<Band>, Box<dyn Error>> {
        self.vectorize(
            VectorizationType::Isoband,
            corners,
            rasterise_scale,
            thresholds,
            interp_method,
            smooth,
        )
    }
}
