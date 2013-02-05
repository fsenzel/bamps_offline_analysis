svn diff svn+ssh://gallmei@th.physik.uni-frankfurt.de/home/bamps/svn/full/offlineAnalysis/trunk/src svn+ssh://gallmei@th.physik.uni-frankfurt.de/home/bamps/svn/full/offlineAnalysis/branches/vector4D/src > diff.txt

svn diff svn+ssh://gallmei@th.physik.uni-frankfurt.de/home/bamps/svn/full/branches/vector4D/src svn+ssh://gallmei@th.physik.uni-frankfurt.de/home/bamps/svn/full/offlineAnalysis/branches/vector4D/src > diff2.txt


'equal' files:

// at revision 902, this file is identical to full/trunk/src/coordinateBins.h
//
// in order to be consistent with the old version, we had to set the
// default of 'timestepScaling = 0.1'(instead of 0.2)

// at revision 902, this file is identical to full/trunk/src/coordinateBins.cpp
//
// every usage of 'particles' is replaced by 'particles_atTimeNow'

// at revision 902, this file is identical to full/branches/vector4D/src/initialmodel.h

// at revision 902 this file is identical to
// full/branches/vector4D/src/initialmodel.cpp

// at revision 1004, this is identical to
// full/branches/vector4D/src/initialmodel_cgc.h

// at revision 902, this file is identical to full/branches/vector4D/src/initialmodel_cgc.cpp

// at revision 902, this file is identical to
full/branches/vector4D/src/offlineoutput.h

// at revision 902, this file is identical to full/branches/vector4D/src/offlineoutput.cpp

// At revision 902, this file is more or less identical to
// full/branches/vector4D/src/particle.h, except N_EVENT_AA and N_EVENT_pp


// at revision 912, this is identical to full/branches/vector4D/src/ringcontainer.h

// at revision 912, this is identical to full/branches/vector4D/src/ringcontainer.cpp


// at revision 902, this file is identical to full/branches/vector4D/src/ringstructure.h
//
// the method 'setLongitudinalGeometry' is new !!!


// at revision 902, this file is identical to full/branches/vector4D/src/ringstructure.cpp
//
// the method 'setLongitudinalGeometry' is new!!!
//
// The original version was written for 'ParticleOffline', but since
// this is derived from 'Particle', we can match it to the latter.

// at revision 1004, this is identical to full/branches/vector4D/src/cellcontainer.h

// at rev 1004, this is identical to
// full/branches/vector4D/src/cellcontainer.cpp,
// except some commenting out in writeAveragesToParticle(...)


funny stuff:

offlineheavyioncollison.cpp:

    int selected = 0;
    do  // pick random particle (that is neihter dead nor the jet) from the _particleList
    {
      selected = static_cast<int>( ran2() * _allParticlesList.size() );
      if ( selected >= _allParticlesList.size() )
      {
        continue;
      }      
      iscat = _allParticlesList[selected];






int coordinateBins::getIndex( const double _x ) const
{
  int index = static_cast<int>(( _x - _min_real ) / _delta_x );

  index -= _min_index_limit;

  if ( index < _min_index_limit )
  {
    index = _min_index_limit;
  }
  if ( index > _max_index_limit )
  {
    index = _max_index_limit;
  }

  return index;
}




          m1 = int ( _allParticlesList.size() * ran2() );
          if ( m1 == _allParticlesList.size() )
          {
            m1 = _allParticlesList.size() - 1;
          }
          iscat = _allParticlesList[m1];
