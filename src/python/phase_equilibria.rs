#[macro_export]
macro_rules! impl_vle_state {
    ($eos:ty, $py_eos:ty) => {
        /// A thermodynamic two phase equilibrium state.
        #[pyclass(name = "PhaseEquilibrium", unsendable)]
        #[derive(Clone)]
        pub struct PyPhaseEquilibrium(PhaseEquilibrium<SIUnit, $eos, 2>);

        #[pymethods]
        impl PyPhaseEquilibrium {
            /// Create a liquid and vapor state in equilibrium
            /// for a pure substance given temperature.
            ///
            /// Parameters
            /// ----------
            /// eos : Saft
            ///     The SAFT equation of state.
            /// temperature : SINumber
            ///     The system temperature.
            /// initial_state : PhaseEquilibrium, optional
            ///     A phase equilibrium used as initial guess.
            ///     Can speed up convergence.
            /// max_iter : int, optional
            ///     The maximum number of iterations.
            /// tol: float, optional
            ///     The solution tolerance.
            /// verbosity : Verbosity, optional
            ///     The verbosity.
            ///
            /// Returns
            /// -------
            /// PhaseEquilibrium
            ///
            /// Raises
            /// ------
            /// RuntimeError
            ///     When pressure iteration fails or no phase equilibrium is found.
            #[staticmethod]
            #[pyo3(text_signature = "(eos, temperature, initial_state=None, max_iter=None, tol=None, verbosity=None)")]
            pub fn pure_t(
                eos: $py_eos,
                temperature: PySINumber,
                initial_state: Option<&PyPhaseEquilibrium>,
                max_iter: Option<usize>,
                tol: Option<f64>,
                verbosity: Option<PyVerbosity>,
            ) -> PyResult<Self> {
                Ok(Self(PhaseEquilibrium::pure_t(
                    &eos.0,
                    temperature.into(),
                    initial_state.and_then(|s| Some(&s.0)),
                    (max_iter, tol, verbosity.map(|v| v.0)).into(),
                )?))
            }

            /// Create a liquid and vapor state in equilibrium
            /// for a pure substance given pressure.
            ///
            /// Parameters
            /// ----------
            /// eos : Saft
            ///     The SAFT equation of state.
            /// pressure : SINumber
            ///     The system pressure.
            /// initial_state : PhaseEquilibrium, optional
            ///     A phase equilibrium used as initial guess.
            ///     Can speed up convergence.
            /// max_iter : int, optional
            ///     The maximum number of iterations.
            /// tol: float, optional
            ///     The solution tolerance.
            /// verbosity : Verbosity, optional
            ///     The verbosity.
            ///
            /// Returns
            /// -------
            /// PhaseEquilibrium
            ///
            /// Raises
            /// ------
            /// RuntimeError
            ///     When pressure iteration fails or no phase equilibrium is found.
            #[staticmethod]
            #[pyo3(text_signature = "(eos, pressure, initial_state=None, max_iter=None, tol=None, verbosity=None)")]
            pub fn pure_p(
                eos: $py_eos,
                pressure: PySINumber,
                initial_state: Option<&PyPhaseEquilibrium>,
                max_iter: Option<usize>,
                tol: Option<f64>,
                verbosity: Option<PyVerbosity>,
            ) -> PyResult<Self> {
                Ok(Self(PhaseEquilibrium::pure_p(
                    &eos.0,
                    pressure.into(),
                    initial_state.and_then(|s| Some(&s.0)),
                    (max_iter, tol, verbosity.map(|v| v.0)).into(),
                )?))
            }

            /// Create a liquid and vapor state in equilibrium
            /// for given temperature, pressure and feed composition.
            ///
            /// Can also be used to calculate liquid liquid phase separation.
            ///
            /// Parameters
            /// ----------
            /// eos : Saft
            ///     The SAFT equation of state.
            /// temperature : SINumber
            ///     The system temperature.
            /// pressure : SINumber
            ///     The system pressure.
            /// feed : SIArray1
            ///     Feed composition (units of amount of substance).
            /// init_vle_state : PhaseEquilibrium, optional
            ///     A phase equilibrium used as initial guess.
            ///     Can speed up convergence.
            /// max_iter : int, optional
            ///     The maximum number of iterations.
            /// tol: float, optional
            ///     The solution tolerance.
            /// verbosity : Verbosity, optional
            ///     The verbosity.
            ///
            /// Returns
            /// -------
            /// PhaseEquilibrium
            ///
            /// Raises
            /// ------
            /// RuntimeError
            ///     When pressure iteration fails or no phase equilibrium is found.
            #[staticmethod]
            #[pyo3(text_signature = "(eos, temperature, pressure, feed, initial_vle_state=None, max_iter=None, tol=None, verbosity=None, non_volatile_components=None)")]
            pub fn tp_flash(
                eos: $py_eos,
                temperature: PySINumber,
                pressure: PySINumber,
                feed: &PySIArray1,
                init_vle_state: Option<&PyPhaseEquilibrium>,
                max_iter: Option<usize>,
                tol: Option<f64>,
                verbosity: Option<PyVerbosity>,
                non_volatile_components: Option<Vec<usize>>,
            ) -> PyResult<Self> {
                Ok(Self(PhaseEquilibrium::tp_flash(
                    &eos.0,
                    temperature.into(),
                    pressure.into(),
                    feed,
                    init_vle_state.and_then(|s| Some(&s.0)),
                    (max_iter, tol, verbosity.map(|v| v.0)).into(), non_volatile_components
                )?))
            }

            /// Compute a VLE given temperature and liquid mole fraction.
            ///
            /// Parameters
            /// ----------
            /// eos : Saft
            ///     The SAFT equation of state.
            /// temperature : SINumber
            ///     The system temperature.
            /// liquid_molefracs : numpy.ndarray
            ///     The mole fraction of the liquid phase.
            /// pressure : SINumber, optional
            ///     The system pressure used as starting condition for iteration.
            /// vapor_molefracs : numpy.ndarray, optional
            ///     The mole fraction of the vapor phase used as
            ///     starting condition for iteration.
            /// max_iter_inner : int, optional
            ///     The maximum number of inner iterations.
            /// max_iter_outer : int, optional
            ///     The maximum number of outer iterations.
            /// tol_inner : float, optional
            ///     The solution tolerance in the inner loop.
            /// tol_outer : float, optional
            ///     The solution tolerance in the outer loop.
            /// verbosity : Verbosity, optional
            ///     The verbosity.
            ///
            /// Returns
            /// -------
            /// PhaseEquilibrium
            #[staticmethod]
            #[pyo3(text_signature = "(eos, temperature, liquid_molefracs, pressure=None, vapor_molefracs=None, max_iter_inner=None, max_iter_outer=None, tol_inner=None, tol_outer=None, verbosity=None)")]
            pub fn bubble_point_tx(
                eos: $py_eos,
                temperature: PySINumber,
                liquid_molefracs: &PyArray1<f64>,
                pressure: Option<PySINumber>,
                vapor_molefracs: Option<&PyArray1<f64>>,
                max_iter_inner: Option<usize>,
                max_iter_outer: Option<usize>,
                tol_inner: Option<f64>,
                tol_outer: Option<f64>,
                verbosity: Option<PyVerbosity>,
            ) -> PyResult<Self> {
                let x = vapor_molefracs.and_then(|m| Some(m.to_owned_array()));
                Ok(Self(PhaseEquilibrium::bubble_point_tx(
                    &eos.0,
                    temperature.into(),
                    pressure.map(|p| p.into()),
                    &liquid_molefracs.to_owned_array(),
                    x.as_ref(),
                    (
                        (max_iter_inner, tol_inner, verbosity.map(|v| v.0)).into(),
                        (max_iter_outer, tol_outer, verbosity.map(|v| v.0)).into()
                    )
                )?))
            }

            /// Compute a VLE given pressure and liquid mole fraction.
            ///
            /// Parameters
            /// ----------
            /// eos : Saft
            ///     The SAFT equation of state.
            /// pressure : SINumber
            ///     The system pressure.
            /// liquid_molefracs : numpy.ndarray
            ///     The mole fraction of the liquid phase.
            /// temperature : SINumber, optional
            ///     The system temperature used as starting condition for iteration.
            /// vapor_molefracs : numpy.ndarray, optional
            ///     The mole fraction of the vapor phase used as
            ///     starting condition for iteration.
            /// max_iter_inner : int, optional
            ///     The maximum number of inner iterations.
            /// max_iter_outer : int, optional
            ///     The maximum number of outer iterations.
            /// tol_inner : float, optional
            ///     The solution tolerance in the inner loop.
            /// tol_outer : float, optional
            ///     The solution tolerance in the outer loop.
            /// verbosity : Verbosity, optional
            ///     The verbosity.
            ///
            /// Returns
            /// -------
            /// PhaseEquilibrium
            #[staticmethod]
            #[pyo3(text_signature = "(eos, pressure, liquid_molefracs, temperature=None, vapor_molefracs=None, max_iter_inner=None, max_iter_outer=None, tol_inner=None, tol_outer=None, verbosity=None)")]
            pub fn bubble_point_px(
                eos: $py_eos,
                pressure: PySINumber,
                liquid_molefracs: &PyArray1<f64>,
                temperature: Option<PySINumber>,
                vapor_molefracs: Option<&PyArray1<f64>>,
                max_iter_inner: Option<usize>,
                max_iter_outer: Option<usize>,
                tol_inner: Option<f64>,
                tol_outer: Option<f64>,
                verbosity: Option<PyVerbosity>,
            ) -> PyResult<Self> {
                let x = vapor_molefracs.and_then(|m| Some(m.to_owned_array()));
                Ok(Self(PhaseEquilibrium::bubble_point_px(
                    &eos.0,
                    pressure.into(),
                    temperature.map(|t| t.into()),
                    &liquid_molefracs.to_owned_array(),
                    x.as_ref(),
                    (
                        (max_iter_inner, tol_inner, verbosity.map(|v| v.0)).into(),
                        (max_iter_outer, tol_outer, verbosity.map(|v| v.0)).into()
                    )
                )?))
            }

            /// Compute a VLE given temperature and vapor mole fraction.
            ///
            /// Parameters
            /// ----------
            /// eos : Saft
            ///     The SAFT equation of state.
            /// temperature : SINumber
            ///     The system temperature.
            /// vapor_molefracs : numpy.ndarray
            ///     The mole fraction of the vapor phase.
            /// pressure : SINumber, optional
            ///     The system pressure used as starting condition for iteration.
            /// liquid_molefracs : numpy.ndarray, optional
            ///     The mole fraction of the liquid phase used as
            ///     starting condition for iteration.
            /// max_iter_inner : int, optional
            ///     The maximum number of inner iterations.
            /// max_iter_outer : int, optional
            ///     The maximum number of outer iterations.
            /// tol_inner : float, optional
            ///     The solution tolerance in the inner loop.
            /// tol_outer : float, optional
            ///     The solution tolerance in the outer loop.
            /// verbosity : Verbosity, optional
            ///     The verbosity.
            ///
            /// Returns
            /// -------
            /// PhaseEquilibrium
            #[staticmethod]
            #[pyo3(text_signature = "(eos, temperature, vapor_molefracs, pressure=None, liquid_molefracs=None, max_iter_inner=None, max_iter_outer=None, tol_inner=None, tol_outer=None, verbosity=None)")]
            pub fn dew_point_tx(
                eos: $py_eos,
                temperature: PySINumber,
                vapor_molefracs: &PyArray1<f64>,
                pressure: Option<PySINumber>,
                liquid_molefracs: Option<&PyArray1<f64>>,
                max_iter_inner: Option<usize>,
                max_iter_outer: Option<usize>,
                tol_inner: Option<f64>,
                tol_outer: Option<f64>,
                verbosity: Option<PyVerbosity>,
            ) -> PyResult<Self> {
                let x = liquid_molefracs.and_then(|m| Some(m.to_owned_array()));
                Ok(Self(PhaseEquilibrium::dew_point_tx(
                    &eos.0,
                    temperature.into(),
                    pressure.map(|p| p.into()),
                    &vapor_molefracs.to_owned_array(),
                    x.as_ref(),
                    (
                        (max_iter_inner, tol_inner, verbosity.map(|v| v.0)).into(),
                        (max_iter_outer, tol_outer, verbosity.map(|v| v.0)).into()
                    )
                )?))
            }

            /// Compute a VLE given pressure and vapor mole fraction.
            ///
            /// Parameters
            /// ----------
            /// eos : Saft
            ///     The SAFT equation of state.
            /// pressure : SINumber
            ///     The system pressure.
            /// liquid_molefracs : numpy.ndarray
            ///     The mole fraction of the liquid phase.
            /// temperature : SINumber, optional
            ///     The system temperature used as starting condition for iteration.
            /// vapor_molefracs : numpy.ndarray, optional
            ///     The mole fraction of the vapor phase used as
            ///     starting condition for iteration.
            /// max_iter_inner : int, optional
            ///     The maximum number of inner iterations.
            /// max_iter_outer : int, optional
            ///     The maximum number of outer iterations.
            /// tol_inner : float, optional
            ///     The solution tolerance in the inner loop.
            /// tol_outer : float, optional
            ///     The solution tolerance in the outer loop.
            /// verbosity : Verbosity, optional
            ///     The verbosity.
            ///
            /// Returns
            /// -------
            /// PhaseEquilibrium
            #[staticmethod]
            #[pyo3(text_signature = "(eos, pressure, vapor_molefracs, temperature=None, liquid_molefracs=None, max_iter_inner=None, max_iter_outer=None, tol_inner=None, tol_outer=None, verbosity=None)")]
            pub fn dew_point_px(
                eos: $py_eos,
                pressure: PySINumber,
                vapor_molefracs: &PyArray1<f64>,
                temperature: Option<PySINumber>,
                liquid_molefracs: Option<&PyArray1<f64>>,
                max_iter_inner: Option<usize>,
                max_iter_outer: Option<usize>,
                tol_inner: Option<f64>,
                tol_outer: Option<f64>,
                verbosity: Option<PyVerbosity>,
            ) -> PyResult<Self> {
                let x = liquid_molefracs.and_then(|m| Some(m.to_owned_array()));
                Ok(Self(PhaseEquilibrium::dew_point_px(
                    &eos.0,
                    pressure.into(),
                    temperature.map(|t| t.into()),
                    &vapor_molefracs.to_owned_array(),
                    x.as_ref(),
                    (
                        (max_iter_inner, tol_inner, verbosity.map(|v| v.0)).into(),
                        (max_iter_outer, tol_outer, verbosity.map(|v| v.0)).into()
                    )
                )?))
            }

            #[getter]
            fn get_vapor(&self) -> PyState {
                PyState(self.0.vapor().clone())
            }

            #[getter]
            fn get_liquid(&self) -> PyState {
                PyState(self.0.liquid().clone())
            }

            /// Calculate a new PhaseEquilibrium with the given chemical potential.
            /// The temperature remains constant, but the states are not in
            /// a mechanical equilibrium anymore.
            ///
            /// Parameters
            /// ----------
            /// chemical_potential: SIArray1
            ///     The new chemical potential
            ///
            #[pyo3(text_signature = "(chemical_potential)")]
            fn update_chemical_potential(slf: &PyCell<Self>, chemical_potential: &PySIArray1) -> PyResult<()> {
                slf.borrow_mut().0.update_chemical_potential(chemical_potential)?;
                Ok(())
            }

            /// Calculate the pure component vapor-liquid equilibria for all
            /// components in the system.
            ///
            /// Parameters
            /// ----------
            /// eos : Saft
            ///     The SAFT equation of state.
            /// temperature : SINumber
            ///     The system temperature.
            ///
            /// Returns
            /// -------
            /// list[PhaseEquilibrium]
            #[staticmethod]
            #[pyo3(text_signature = "(eos, temperature)")]
            fn vle_pure_comps_t(eos: $py_eos, temperature: PySINumber) -> Vec<Option<Self>> {
                PhaseEquilibrium::vle_pure_comps_t(&eos.0, temperature.into())
                    .into_iter()
                    .map(|o| o.map(Self))
                    .collect()
            }

            /// Calculate the pure component vapor-liquid equilibria for all
            /// components in the system.
            ///
            /// Parameters
            /// ----------
            /// eos : Saft
            ///     The SAFT equation of state.
            /// pressure : SINumber
            ///     The system pressure.
            ///
            /// Returns
            /// -------
            /// list[PhaseEquilibrium]
            #[staticmethod]
            #[pyo3(text_signature = "(eos, pressure)")]
            fn vle_pure_comps_p(eos: $py_eos, pressure: PySINumber) -> Vec<Option<Self>> {
                PhaseEquilibrium::vle_pure_comps_p(&eos.0, pressure.into())
                    .into_iter()
                    .map(|o| o.map(Self))
                    .collect()
            }

            /// Calculate the pure component vapor pressures for all the
            /// components in the system.
            ///
            /// Parameters
            /// ----------
            /// eos : Saft
            ///     The SAFT equation of state.
            /// temperature : SINumber
            ///     The system temperature.
            ///
            /// Returns
            /// -------
            /// list[SINumber]
            #[staticmethod]
            #[pyo3(text_signature = "(eos, temperature)")]
            fn vapor_pressure(eos: $py_eos, temperature: PySINumber) -> Vec<Option<PySINumber>> {
                PhaseEquilibrium::vapor_pressure(&eos.0, temperature.into())
                    .into_iter()
                    .map(|o| o.map(|n| n.into()))
                    .collect()
            }

            /// Calculate the pure component boiling temperatures for all the
            /// components in the system.
            ///
            /// Parameters
            /// ----------
            /// eos : Saft
            ///     The SAFT equation of state.
            /// pressure : SINumber
            ///     The system pressure.
            ///
            /// Returns
            /// -------
            /// list[SINumber]
            #[staticmethod]
            #[pyo3(text_signature = "(eos, pressure)")]
            fn boiling_temperature(eos: $py_eos, pressure: PySINumber) -> Vec<Option<PySINumber>> {
                PhaseEquilibrium::boiling_temperature(&eos.0, pressure.into())
                    .into_iter()
                    .map(|o| o.map(|n| n.into()))
                    .collect()
            }

            fn _repr_markdown_(&self) -> String {
                self.0._repr_markdown_()
            }
        }

        #[pyproto]
        impl pyo3::class::basic::PyObjectProtocol for PyPhaseEquilibrium {
            fn __repr__(&self) -> PyResult<String> {
                Ok(self.0.to_string())
            }
        }

        /// A thermodynamic three phase equilibrium state.
        #[pyclass(name = "ThreePhaseEquilibrium", unsendable)]
        #[derive(Clone)]
        struct PyThreePhaseEquilibrium(PhaseEquilibrium<SIUnit, $eos, 3>);

        #[pymethods]
        impl PyPhaseEquilibrium {
            /// Calculate a heteroazeotrope in a binary mixture for a given temperature.
            ///
            /// Parameters
            /// ----------
            /// eos : Saft
            ///     The SAFT equation of state.
            /// temperature : SINumber
            ///     The system temperature.
            /// x_init : list[float]
            ///     Initial guesses for the liquid molefracs of component 1
            ///     at the heteroazeotropic point.
            /// max_iter : int, optional
            ///     The maximum number of iterations.
            /// tol: float, optional
            ///     The solution tolerance.
            /// verbosity : Verbosity, optional
            ///     The verbosity.
            /// max_iter_bd_inner : int, optional
            ///     The maximum number of inner iterations in the bubble/dew point iteration.
            /// max_iter_bd_outer : int, optional
            ///     The maximum number of outer iterations in the bubble/dew point iteration.
            /// tol_bd_inner : float, optional
            ///     The solution tolerance in the inner loop of the bubble/dew point iteration.
            /// tol_bd_outer : float, optional
            ///     The solution tolerance in the outer loop of the bubble/dew point iteration.
            /// verbosity_bd : Verbosity, optional
            ///     The verbosity of the bubble/dew point iteration.
            #[staticmethod]
            #[pyo3(text_signature = "(eos, temperature, x_init, max_iter=None, tol=None, verbosity=None, max_iter_bd_inner=None, max_iter_bd_outer=None, tol_bd_inner=None, tol_bd_outer=None, verbosity_bd=None)")]
            fn heteroazeotrope_t(
                eos: $py_eos,
                temperature: PySINumber,
                x_init: (f64, f64),
                max_iter: Option<usize>,
                tol: Option<f64>,
                verbosity: Option<PyVerbosity>,
                max_iter_bd_inner: Option<usize>,
                max_iter_bd_outer: Option<usize>,
                tol_bd_inner: Option<f64>,
                tol_bd_outer: Option<f64>,
                verbosity_bd: Option<PyVerbosity>,
            ) -> PyResult<PyThreePhaseEquilibrium> {
                Ok(PyThreePhaseEquilibrium(PhaseEquilibrium::heteroazeotrope_t(
                    &eos.0,
                    temperature.into(),
                    x_init,
                    (max_iter, tol, verbosity.map(|v| v.0)).into(),
                    (
                        (max_iter_bd_inner, tol_bd_inner, verbosity_bd.map(|v| v.0)).into(),
                        (max_iter_bd_outer, tol_bd_outer, verbosity_bd.map(|v| v.0)).into(),
                    )
                )?))
            }

            /// Calculate a heteroazeotrope in a binary mixture for a given pressure.
            ///
            /// Parameters
            /// ----------
            /// eos : Saft
            ///     The SAFT equation of state.
            /// pressure : SINumber
            ///     The system pressure.
            /// x_init : list[float]
            ///     Initial guesses for the liquid molefracs of component 1
            ///     at the heteroazeotropic point.
            /// max_iter : int, optional
            ///     The maximum number of iterations.
            /// tol: float, optional
            ///     The solution tolerance.
            /// verbosity : Verbosity, optional
            ///     The verbosity.
            /// max_iter_bd_inner : int, optional
            ///     The maximum number of inner iterations in the bubble/dew point iteration.
            /// max_iter_bd_outer : int, optional
            ///     The maximum number of outer iterations in the bubble/dew point iteration.
            /// tol_bd_inner : float, optional
            ///     The solution tolerance in the inner loop of the bubble/dew point iteration.
            /// tol_bd_outer : float, optional
            ///     The solution tolerance in the outer loop of the bubble/dew point iteration.
            /// verbosity_bd : Verbosity, optional
            ///     The verbosity of the bubble/dew point iteration.
            #[staticmethod]
            #[pyo3(text_signature = "(eos, pressure, x_init, max_iter=None, tol=None, verbosity=None, max_iter_bd_inner=None, max_iter_bd_outer=None, tol_bd_inner=None, tol_bd_outer=None, verbosity_bd=None)")]
            fn heteroazeotrope_p(
                eos: $py_eos,
                pressure: PySINumber,
                x_init: (f64, f64),
                max_iter: Option<usize>,
                tol: Option<f64>,
                verbosity: Option<PyVerbosity>,
                max_iter_bd_inner: Option<usize>,
                max_iter_bd_outer: Option<usize>,
                tol_bd_inner: Option<f64>,
                tol_bd_outer: Option<f64>,
                verbosity_bd: Option<PyVerbosity>,
            ) -> PyResult<PyThreePhaseEquilibrium> {
                Ok(PyThreePhaseEquilibrium(PhaseEquilibrium::heteroazeotrope_p(
                    &eos.0,
                    pressure.into(),
                    x_init,
                    (max_iter, tol, verbosity.map(|v| v.0)).into(),
                    (
                        (max_iter_bd_inner, tol_bd_inner, verbosity_bd.map(|v| v.0)).into(),
                        (max_iter_bd_outer, tol_bd_outer, verbosity_bd.map(|v| v.0)).into(),
                    )
                )?))
            }
        }

        #[pymethods]
        impl PyThreePhaseEquilibrium {
            #[getter]
            fn get_vapor(&self) -> PyState {
                PyState(self.0.vapor().clone())
            }

            #[getter]
            fn get_liquid1(&self) -> PyState {
                PyState(self.0.liquid1().clone())
            }

            #[getter]
            fn get_liquid2(&self) -> PyState {
                PyState(self.0.liquid2().clone())
            }

            fn _repr_markdown_(&self) -> String {
                self.0._repr_markdown_()
            }
        }

        #[pyproto]
        impl pyo3::class::basic::PyObjectProtocol for PyThreePhaseEquilibrium {
            fn __repr__(&self) -> PyResult<String> {
                Ok(self.0.to_string())
            }
        }

        #[pymethods]
        impl PyState {
            /// Calculates a two phase Tp-flash with the state as feed.
            ///
            /// Parameters
            /// ----------
            /// init_vle_state : PhaseEquilibrium, optional
            ///     A phase equilibrium used as initial guess.
            ///     Can speed up convergence.
            /// max_iter : int, optional
            ///     The maximum number of iterations.
            /// tol: float, optional
            ///     The solution tolerance.
            /// verbosity : Verbosity, optional
            ///     The verbosity.
            ///
            /// Returns
            /// -------
            /// PhaseEquilibrium
            ///
            /// Raises
            /// ------
            /// RuntimeError
            ///     When pressure iteration fails or no phase equilibrium is found.
            #[pyo3(text_signature = "($self, initial_vle_state=None, max_iter=None, tol=None, verbosity=None, non_volatile_components=None)")]
            pub fn tp_flash(
                &self,
                init_vle_state: Option<&PyPhaseEquilibrium>,
                max_iter: Option<usize>,
                tol: Option<f64>,
                verbosity: Option<PyVerbosity>,
                non_volatile_components: Option<Vec<usize>>,
            ) -> PyResult<PyPhaseEquilibrium> {
                Ok(PyPhaseEquilibrium(self.0.tp_flash(
                    init_vle_state.and_then(|s| Some(&s.0)),
                    (max_iter, tol, verbosity.map(|v| v.0)).into(),
                    non_volatile_components
                )?))
            }
        }

        /// Phase diagram for a pure component.
        ///
        /// Parameters
        /// ----------
        /// eos: Eos
        ///     The equation of state.
        /// min_temperature: SINumber
        ///     The lower limit for the temperature.
        /// npoints: int
        ///     The number of points.
        /// critical_temperature: SINumber, optional
        ///     An estimate for the critical temperature to initialize
        ///     the calculation if necessary. For most components not necessary.
        ///     Defaults to `None`.
        /// max_iter : int, optional
        ///     The maximum number of iterations.
        /// tol: float, optional
        ///     The solution tolerance.
        /// verbosity : Verbosity, optional
        ///     The verbosity.
        ///
        /// Returns
        /// -------
        /// PhaseDiagramPure
        #[pyclass(name = "PhaseDiagramPure", unsendable)]
        #[pyo3(text_signature = "(eos, min_temperature, npoints, critical_temperature=None, max_iter=None, tol=None, verbosity=None)")]
        pub struct PyPhaseDiagramPure(PhaseDiagramPure<SIUnit, $eos>);

        #[pymethods]
        impl PyPhaseDiagramPure {
            #[new]
            pub fn new(
                eos: &$py_eos,
                min_temperature: PySINumber,
                npoints: usize,
                critical_temperature: Option<PySINumber>,
                max_iter: Option<usize>,
                tol: Option<f64>,
                verbosity: Option<PyVerbosity>,
            ) -> PyResult<Self> {
                let dia = PhaseDiagramPure::new(
                    &eos.0,
                    min_temperature.into(),
                    npoints,
                    critical_temperature.map(|t| t.into()),
                    (max_iter, tol, verbosity.map(|v| v.0)).into(),
                )?;
                Ok(Self(dia))
            }

            #[getter]
            pub fn get_states(&self) -> Vec<PyPhaseEquilibrium> {
                self.0
                    .states
                    .iter()
                    .map(|vle| PyPhaseEquilibrium(vle.clone()))
                    .collect()
            }

            #[getter]
            pub fn get_temperature(&self) -> PySIArray1 {
                self.0.temperature().into()
            }

            #[getter]
            pub fn get_pressure(&self) -> PySIArray1 {
                self.0.pressure().into()
            }

            #[getter]
            pub fn get_density_vapor(&self) -> PySIArray1 {
                self.0.density_vapor().into()
            }

            #[getter]
            pub fn get_density_liquid(&self) -> PySIArray1 {
                self.0.density_liquid().into()
            }

            #[getter]
            pub fn get_molar_enthalpy_vapor(&self) -> PySIArray1 {
                self.0.molar_enthalpy_vapor().into()
            }

            #[getter]
            pub fn get_molar_enthalpy_liquid(&self) -> PySIArray1 {
                self.0.molar_enthalpy_liquid().into()
            }

            #[getter]
            pub fn get_molar_entropy_vapor(&self) -> PySIArray1 {
                self.0.molar_entropy_vapor().into()
            }

            #[getter]
            pub fn get_molar_entropy_liquid(&self) -> PySIArray1 {
                self.0.molar_entropy_liquid().into()
            }

            /// Returns the phase diagram as dictionary.
            ///
            /// Note
            /// ----
            /// temperature : K
            /// pressure : Pa
            /// densities : mol / m³
            /// molar enthalpies : kJ / mol
            /// molar entropies : kJ / mol / K
            ///
            /// Returns
            /// -------
            /// dict[str, list[float]]
            ///     Keys: property names. Values: property for each state.
            pub fn to_dict(&self) -> PyResult<HashMap<String, Vec<f64>>> {
                let mut result = HashMap::with_capacity(8);
                result.insert(String::from("temperature"), (self.0.temperature() / KELVIN).into_value()?.into_raw_vec());
                result.insert(String::from("pressure"), (self.0.pressure() / PASCAL).into_value()?.into_raw_vec());
                result.insert(String::from("density liquid"), (self.0.density_liquid() / (MOL / METER.powi(3))).into_value()?.into_raw_vec());
                result.insert(String::from("density vapor"), (self.0.density_vapor() / (MOL / METER.powi(3))).into_value()?.into_raw_vec());
                result.insert(String::from("molar enthalpy liquid"), (self.0.molar_enthalpy_liquid() / (KILO*JOULE / MOL)).into_value()?.into_raw_vec());
                result.insert(String::from("molar enthalpy vapor"), (self.0.molar_enthalpy_vapor() / (KILO*JOULE / MOL)).into_value()?.into_raw_vec());
                result.insert(String::from("molar entropy liquid"), (self.0.molar_entropy_liquid() / (KILO*JOULE / KELVIN / MOL)).into_value()?.into_raw_vec());
                result.insert(String::from("molar entropy vapor"), (self.0.molar_entropy_vapor() / (KILO*JOULE / KELVIN / MOL)).into_value()?.into_raw_vec());
                Ok(result)
            }

        }


        /// Phase diagram for a binary mixture.
        #[pyclass(name = "PhaseDiagramBinary", unsendable)]
        pub struct PyPhaseDiagramBinary(PhaseDiagramBinary<SIUnit, $eos>);

        #[pymethods]
        impl PyPhaseDiagramBinary {
            /// Txy phase diagram for a binary mixture.
            ///
            /// Parameters
            /// ----------
            /// eos: SaftFunctional
            ///     The SAFT Helmholtz energy functional.
            /// pressure: SINumber
            ///     The pressure.
            /// npoints: int, optional
            ///     The number of points (default 51).
            /// x_lle: SINumber, optional
            ///     An estimate for the molefractions of component 1
            ///     at the heteroazeotrop
            /// max_iter_inner : int, optional
            ///     The maximum number of inner iterations in the bubble/dew point iteration.
            /// max_iter_outer : int, optional
            ///     The maximum number of outer iterations in the bubble/dew point iteration.
            /// tol_inner : float, optional
            ///     The solution tolerance in the inner loop of the bubble/dew point iteration.
            /// tol_outer : float, optional
            ///     The solution tolerance in the outer loop of the bubble/dew point iteration.
            /// verbosity : Verbosity, optional
            ///     The verbosity of the bubble/dew point iteration.
            ///
            /// Returns
            /// -------
            /// PhaseDiagramBinary
            #[staticmethod]
            #[pyo3(text_signature = "(eos, pressure, npoints=None, x_lle=None, max_iter_bd_inner=None, max_iter_bd_outer=None, tol_bd_inner=None, tol_bd_outer=None, verbosity_bd=None)")]
            pub fn new_txy(
                eos: $py_eos,
                pressure: PySINumber,
                npoints: Option<usize>,
                x_lle: Option<(f64, f64)>,
                max_iter_inner: Option<usize>,
                max_iter_outer: Option<usize>,
                tol_inner: Option<f64>,
                tol_outer: Option<f64>,
                verbosity: Option<PyVerbosity>,
            ) -> PyResult<Self> {
                let dia = PhaseDiagramBinary::new_txy(
                    &eos.0,
                    pressure.into(),
                    npoints,
                    x_lle,
                    (
                        (max_iter_inner, tol_inner, verbosity.map(|v| v.0)).into(),
                        (max_iter_outer, tol_outer, verbosity.map(|v| v.0)).into(),
                    )
                )?;
                Ok(Self(dia))
            }

            /// pxy phase diagram for a binary mixture.
            ///
            /// Parameters
            /// ----------
            /// eos: SaftFunctional
            ///     The SAFT Helmholtz energy functional.
            /// temperature: SINumber
            ///     The temperature.
            /// npoints: int, optional
            ///     The number of points (default 51).
            /// x_lle: SINumber, optional
            ///     An estimate for the molefractions of component 1
            ///     at the heteroazeotrop
            /// max_iter_inner : int, optional
            ///     The maximum number of inner iterations in the bubble/dew point iteration.
            /// max_iter_outer : int, optional
            ///     The maximum number of outer iterations in the bubble/dew point iteration.
            /// tol_inner : float, optional
            ///     The solution tolerance in the inner loop of the bubble/dew point iteration.
            /// tol_outer : float, optional
            ///     The solution tolerance in the outer loop of the bubble/dew point iteration.
            /// verbosity : Verbosity, optional
            ///     The verbosity of the bubble/dew point iteration.
            ///
            /// Returns
            /// -------
            /// PhaseDiagramBinary
            #[staticmethod]
            #[pyo3(text_signature = "(eos, temperature, npoints=None, x_lle=None, max_iter_bd_inner=None, max_iter_bd_outer=None, tol_bd_inner=None, tol_bd_outer=None, verbosity_bd=None)")]
            pub fn new_pxy(
                eos: $py_eos,
                temperature: PySINumber,
                npoints: Option<usize>,
                x_lle: Option<(f64, f64)>,
                max_iter_inner: Option<usize>,
                max_iter_outer: Option<usize>,
                tol_inner: Option<f64>,
                tol_outer: Option<f64>,
                verbosity: Option<PyVerbosity>,
            ) -> PyResult<Self> {
                let dia = PhaseDiagramBinary::new_pxy(
                    &eos.0,
                    temperature.into(),
                    npoints,
                    x_lle,
                    (
                        (max_iter_inner, tol_inner, verbosity.map(|v| v.0)).into(),
                        (max_iter_outer, tol_outer, verbosity.map(|v| v.0)).into(),
                    )
                )?;
                Ok(Self(dia))
            }

            /// Txy phase diagram for a liquid-liquid equilibrium of a binary mixture.
            ///
            /// Parameters
            /// ----------
            /// eos: SaftFunctional
            ///     The SAFT Helmholtz energy functional.
            /// pressure: SINumber
            ///     The pressure.
            /// x_feed: float
            ///     Molefraction of component 1 in the (unstable) feed state.
            /// min_temperature:
            ///     The lower limit of the temperature range.
            /// max_temperature:
            ///     The upper limit of the temperature range.
            /// npoints: int, optional
            ///     The number of points (default 51).
            ///
            /// Returns
            /// -------
            /// PhaseDiagramBinary
            #[staticmethod]
            #[pyo3(text_signature = "(eos, pressure, min_temperature, max_temperature, npoints=None)")]
            pub fn new_txy_lle(
                eos: $py_eos,
                pressure: PySINumber,
                x_feed: f64,
                min_temperature: PySINumber,
                max_temperature: PySINumber,
                npoints: Option<usize>,
            ) -> PyResult<Self> {
                let dia = PhaseDiagramBinary::new_txy_lle(
                    &eos.0,
                    pressure.into(),
                    x_feed,
                    min_temperature.into(),
                    max_temperature.into(),
                    npoints,
                )?;
                Ok(Self(dia))
            }

            /// pxy phase diagram for a liquid-liquid equilibrium of a binary mixture.
            ///
            /// Parameters
            /// ----------
            /// eos: SaftFunctional
            ///     The SAFT Helmholtz energy functional.
            /// temperature: SINumber
            ///     The temperature.
            /// x_feed: float
            ///     Molefraction of component 1 in the (unstable) feed state.
            /// min_pressure:
            ///     The lower limit of the pressure range.
            /// max_pressure:
            ///     The upper limit of the pressure range.
            /// npoints: int, optional
            ///     The number of points (default 51).
            ///
            /// Returns
            /// -------
            /// PhaseDiagramBinary
            #[staticmethod]
            #[pyo3(text_signature = "(eos, temperature, min_pressure, max_pressure, npoints=None)")]
            pub fn new_pxy_lle(
                eos: $py_eos,
                temperature: PySINumber,
                x_feed: f64,
                min_pressure: PySINumber,
                max_pressure: PySINumber,
                npoints: Option<usize>,
            ) -> PyResult<Self> {
                let dia = PhaseDiagramBinary::new_pxy_lle(
                    &eos.0,
                    temperature.into(),
                    x_feed,
                    min_pressure.into(),
                    max_pressure.into(),
                    npoints,
                )?;
                Ok(Self(dia))
            }

            #[getter]
            pub fn get_states(&self) -> Vec<PyPhaseEquilibrium> {
                self.0
                    .states
                    .iter()
                    .map(|vle| PyPhaseEquilibrium(vle.clone()))
                    .collect()
            }

            #[getter]
            pub fn get_temperature(&self) -> PySIArray1 {
                self.0.temperature().into()
            }

            #[getter]
            pub fn get_pressure(&self) -> PySIArray1 {
                self.0.pressure().into()
            }

            #[getter]
            fn get_vapor_molefracs<'py>(&self, py: Python<'py>) -> &'py PyArray1<f64> {
                self.0.vapor_molefracs().view().to_pyarray(py)
            }

            #[getter]
            fn get_liquid_molefracs<'py>(&self, py: Python<'py>) -> &'py PyArray1<f64> {
                self.0.liquid_molefracs().view().to_pyarray(py)
            }
        }

        /// Phase diagram for a binary mixture exhibiting a heteroazeotrope.
        #[pyclass(name = "PhaseDiagramHetero", unsendable)]
        pub struct PyPhaseDiagramHetero(PhaseDiagramHetero<SIUnit, $eos>);

        #[pymethods]
        impl PyPhaseDiagramHetero {
            /// Txy phase diagram for a binary mixture exhibiting a heteroazeotrope.
            ///
            /// Parameters
            /// ----------
            /// eos: SaftFunctional
            ///     The SAFT Helmholtz energy functional.
            /// pressure: SINumber
            ///     The pressure.
            /// x_lle: SINumber
            ///     Initial values for the molefractions of component 1
            ///     at the heteroazeotrop.
            /// min_temperature_lle: SINumber, optional
            ///     The minimum temperature up to which the LLE is calculated.
            ///     If it is not provided, no LLE is calcualted.
            /// npoints_vle: int, optional
            ///     The number of points for the VLE (default 51).
            /// npoints_lle: int, optional
            ///     The number of points for the LLE (default 51).
            /// max_iter_inner : int, optional
            ///     The maximum number of inner iterations in the bubble/dew point iteration.
            /// max_iter_outer : int, optional
            ///     The maximum number of outer iterations in the bubble/dew point iteration.
            /// tol_inner : float, optional
            ///     The solution tolerance in the inner loop of the bubble/dew point iteration.
            /// tol_outer : float, optional
            ///     The solution tolerance in the outer loop of the bubble/dew point iteration.
            /// verbosity : Verbosity, optional
            ///     The verbosity of the bubble/dew point iteration.
            ///
            /// Returns
            /// -------
            /// PhaseDiagramBinary
            #[staticmethod]
            #[pyo3(text_signature = "(eos, pressure, x_lle, min_temperature_lle=None, npoints_vle=None, npoints_lle=None, max_iter_bd_inner=None, max_iter_bd_outer=None, tol_bd_inner=None, tol_bd_outer=None, verbosity_bd=None)")]
            pub fn new_txy(
                eos: $py_eos,
                pressure: PySINumber,
                x_lle: (f64, f64),
                min_temperature_lle: Option<PySINumber>,
                npoints_vle: Option<usize>,
                npoints_lle: Option<usize>,
                max_iter_inner: Option<usize>,
                max_iter_outer: Option<usize>,
                tol_inner: Option<f64>,
                tol_outer: Option<f64>,
                verbosity: Option<PyVerbosity>,
            ) -> PyResult<Self> {
                let dia = PhaseDiagramHetero::new_txy(
                    &eos.0,
                    pressure.into(),
                    x_lle,
                    min_temperature_lle.map(|t| t.into()),
                    npoints_vle,
                    npoints_lle,
                    (
                        (max_iter_inner, tol_inner, verbosity.map(|v| v.0)).into(),
                        (max_iter_outer, tol_outer, verbosity.map(|v| v.0)).into(),
                    )
                )?;
                Ok(Self(dia))
            }

            /// pxy phase diagram for a binary mixture exhibiting a heteroazeotrope.
            ///
            /// Parameters
            /// ----------
            /// eos: SaftFunctional
            ///     The SAFT Helmholtz energy functional.
            /// temperature: SINumber
            ///     The temperature.
            /// x_lle: SINumber
            ///     Initial values for the molefractions of component 1
            ///     at the heteroazeotrop.
            /// max_pressure: SINumber, optional
            ///     The maximum pressure up to which the LLE is calculated.
            ///     If it is not provided, no LLE is calcualted.
            /// npoints_vle: int, optional
            ///     The number of points for the VLE (default 51).
            /// npoints_lle: int, optional
            ///     The number of points for the LLE (default 51).
            /// max_iter_inner : int, optional
            ///     The maximum number of inner iterations in the bubble/dew point iteration.
            /// max_iter_outer : int, optional
            ///     The maximum number of outer iterations in the bubble/dew point iteration.
            /// tol_inner : float, optional
            ///     The solution tolerance in the inner loop of the bubble/dew point iteration.
            /// tol_outer : float, optional
            ///     The solution tolerance in the outer loop of the bubble/dew point iteration.
            /// verbosity : Verbosity, optional
            ///     The verbosity of the bubble/dew point iteration.
            ///
            /// Returns
            /// -------
            /// PhaseDiagramBinary
            #[staticmethod]
            #[pyo3(text_signature = "(eos, temperature, x_lle, max_pressure_lle=None, npoints_vle=None, npoints_lle=None, max_iter_bd_inner=None, max_iter_bd_outer=None, tol_bd_inner=None, tol_bd_outer=None, verbosity_bd=None)")]
            pub fn new_pxy(
                eos: $py_eos,
                temperature: PySINumber,
                x_lle: (f64, f64),
                max_pressure_lle: Option<PySINumber>,
                npoints_vle: Option<usize>,
                npoints_lle: Option<usize>,
                max_iter_inner: Option<usize>,
                max_iter_outer: Option<usize>,
                tol_inner: Option<f64>,
                tol_outer: Option<f64>,
                verbosity: Option<PyVerbosity>,
            ) -> PyResult<Self> {
                let dia = PhaseDiagramHetero::new_pxy(
                    &eos.0,
                    temperature.into(),
                    x_lle,
                    max_pressure_lle.map(|t| t.into()),
                    npoints_vle,
                    npoints_lle,
                    (
                        (max_iter_inner, tol_inner, verbosity.map(|v| v.0)).into(),
                        (max_iter_outer, tol_outer, verbosity.map(|v| v.0)).into(),
                    )
                )?;
                Ok(Self(dia))
            }

            #[getter]
            pub fn get_vle(&self) -> PyPhaseDiagramBinary {
                PyPhaseDiagramBinary(self.0.vle().clone())
            }

            #[getter]
            pub fn get_vle1(&self) -> PyPhaseDiagramBinary {
                PyPhaseDiagramBinary(self.0.vle1.clone())
            }

            #[getter]
            pub fn get_vle2(&self) -> PyPhaseDiagramBinary {
                PyPhaseDiagramBinary(self.0.vle2.clone())
            }

            #[getter]
            pub fn get_lle(&self) -> Option<PyPhaseDiagramBinary> {
                self.0
                    .lle
                    .as_ref()
                    .map(|d| PyPhaseDiagramBinary(d.clone()))
            }
        }
    }
}
