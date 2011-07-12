SET(project_SRCS
additionalparticlesdistribution.cpp
analysis.cpp
binary_cross_sections.cpp
binning.cpp
binning2.cpp
cellcontainer.cpp
configuration.cpp
coordinateBins.cpp
FPT_compare.cpp
initialstatemodel.cpp
integrand_time.cpp
integrand23.cpp
integrand32.cpp
integrate.cpp
interactiontype.cpp
interpolation.cpp
interpolation23.cpp
# interpolation_gQ.cpp
lorentz.cpp
minijets.cpp
offlineBAMPS.cpp
offlineheavyioncollison.cpp 
particle.cpp
prefactors23.cpp
psf.cpp
pythiaInitialDistribution.cpp
random.cpp
ratesmanager.cpp
ringcontainer.cpp
ringstructure.cpp
scattering22.cpp
scattering23.cpp
scattering32.cpp
setting.cpp
vegas.cpp
woodsaxon.cpp
)

ADD_EXECUTABLE(offlineBAMPS ${project_SRCS})
#TARGET_LINK_LIBRARIES(offlineBAMPS ${Boost_LIBRARIES})
#TARGET_LINK_LIBRARIES(offlineBAMPS ${Boost_LIBRARIES} ${CUBA_LIBRARY})

#install(TARGETS offlineBAMPS ${INSTALL_TARGETS_DEFAULT_ARGS})
