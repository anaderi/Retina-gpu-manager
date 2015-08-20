#include <cmath>
#include <iostream>
#include <algorithm>
#include <map>
#include <set>
#include <cstdlib>

#include "Retina.h"
#include "Grid.h"
#include "Physics.h"
#include "Tools.h"
#include "Logger.h"

TrackPure generateTrackFromIndex(
  const std::vector<int>& multiIndex,
  const std::vector<std::vector<double> > & dim
)
{
  return TrackPure(
    dim[0][multiIndex[0]],
    dim[1][multiIndex[1]],
    dim[2][multiIndex[2]],
    dim[3][multiIndex[3]]
  );
}

std::vector<TrackPure> restoreTracks(
  const std::vector<std::vector<double> >& dimensions,
  const std::vector<TrackPure>& grid,
  const std::vector<double>& responces
)
{
  int gridSizeNoBoarders = calculateGridSize(dimensions, 1);
  std::vector<int> indexes(dimensions.size(), 1);
  std::vector<TrackPure> restored;
  
  const double threshold = 1.9;//getQuatile(responces, TRACK_THRESHOLD);
  for (int i = 0; i < gridSizeNoBoarders; ++i)
  {
    std::vector<int> neighbours = generateNeighboursIndexes(indexes, dimensions);
    int currentIndex = multiIndexToIndex(indexes, dimensions);
    bool isLocalMaximum = true;
    for (int neighbour : neighbours)
    {
      if (responces[neighbour] > responces[currentIndex])
      {
        isLocalMaximum = false;
      }
    }
    isLocalMaximum = /*isLocalMaxumum ||*/ (responces[currentIndex] > threshold);
    if (isLocalMaximum && responces[currentIndex] > 1e-5)
    {
      TrackPure answer = grid[currentIndex] * responces[currentIndex];
      double sum_responce = responces[currentIndex];
      for (int neighbour : neighbours)
      {
        answer = answer + grid[neighbour] * responces[neighbour];
        sum_responce += responces[neighbour];
      }
      //std::cerr << responces[currentIndex] << "->" << sum_responce << std::endl;
      const TrackPure track = answer * (1.0 / sum_responce);
      restored.push_back(track);      
    /*
     *
     * 
     *  std::cerr << track.xOnZ0 << "," << track.yOnZ0 
        << "," << track.dxOverDz << "," 
        << track.dyOverDz << std::endl; 
     */
    }
    generateNextMultiIndex(indexes, dimensions, 1);
  }
  return restored;
}

std::vector<double> calculateResponses(
  const std::vector<TrackPure>& grid,
  const std::vector<Hit>& hits,
  double sharpness
)
{
  std::vector<double> responces(grid.size());
  for (unsigned int i = 0; i < grid.size(); ++i) 
  {
    auto& track = grid[i];
    for (const Hit& hit : hits)
    {
      responces[i] += exp(-getDistance(track, hit) / sharpness);
    }
  }
  return responces;
}

std::vector<Track> findHits(
  const std::vector<TrackPure>& tracks,
  const std::vector<Hit>& hits
)
{
  std::vector<Track> extendedTracks;
  extendedTracks.reserve(tracks.size());
  std::set<std::vector<uint32_t> > tracksSet;
  std::vector<bool> used(hits.size(), false);
  for (const TrackPure& track: tracks)
  {

    std::map<uint32_t, Hit> sensorsBest;
    for (size_t i = 0; i < hits.size(); ++i)
    {
      //if (!used[i])
      {
        const Hit& hit = hits[i];
        if (!sensorsBest.count(hit.sensorId) ||
            (getDistance(track, sensorsBest[hit.sensorId]) > 
            getDistance(track, hit)))
        {
          sensorsBest[hit.sensorId] = hit;
        }
      }
    }
    std::vector<std::pair<double, Hit> > distances;
    Track extended;
    for (const auto& pair: sensorsBest)
    {
      distances.emplace_back(getDistance(track, pair.second), pair.second);
    }
    std::sort(distances.begin(), distances.end(), 
      [](const std::pair<double, Hit>& a, const std::pair<double, Hit>& b) -> bool
      {
        return a.first < b.first;
      }
    );
    {
      double zStart = distances[0].second.z;
      extended.addHit(distances[0].second.id);
      for (size_t i = 1; i < 3; /*distances.size()*/ i++)
      {
        //if (isFit(track, distances[i].second, zStart))
        {
          extended.addHit(distances[i].second.id);
        }
      }
    }
    if (extended.hitsNum > 2)
    {
      /*for (size_t i = 0; i < extended.hitsNum; ++i)
      {
        used[extended.hits[i]] = true;
      }*/
      extendedTracks.push_back(extended);
    }
  }
        
  return extendedTracks;
}
TrackPure makeLocalMaximum(TrackPure track, const std::vector<Hit>& hits, double rs) {
  std::vector<TrackPure> steps = {
    TrackPure(1, 0, 0, 0),
    TrackPure(0, 1, 0, 0),
    TrackPure(0, 0, 1, 0),
    TrackPure(0, 0, 0, 1)
  };
  double len = 1e-3;
  while (len > 1e-7)
  {
    bool update = false;
    TrackPure best = track;
    double bestResponce = calculateResponses({track}, hits, rs)[0];
    for (double dx: {-len, len})
    {
      for (TrackPure step: steps)
      {
        double current = calculateResponses({track + step * dx}, hits, rs)[0];
        if ( current > bestResponce + 1e-7)
        {
          update = true;
          best = track + step * dx;
          bestResponce = current;
        }
      }
    }
    if (!update)
    {
      len /= 2;
    }
    else
    {
      track = best;
      len *= 2;
    }
    //std::cerr << len << std::endl;
    //std::cerr << bestResponce << std::endl;
  }
  return track;
}
/**
 * Common entrypoint for Gaudi and non-Gaudi
 * @param input  
 * @param output 
 */
int cpuRetinaInvocation(
    const std::vector<const std::vector<uint8_t>* > & input,
    std::vector<std::vector<uint8_t> > & output
)
{
  output.resize(input.size());
  for (size_t i = 0; i < input.size(); ++i)
  {
    EventInfo event = parseEvent(const_cast<const uint8_t*>(&(*input[i])[0]), input[i]->size());
    const std::vector<std::vector<double> > dimensions = generateDimensions(event);
    const std::vector<TrackPure> grid = generateGrid<TrackPure>(dimensions, generateTrackFromIndex);
    const std::vector<Hit>& hits = event.hits;
    //const std::vector<double> responses = calculateResponses(grid, hits, RETINA_SHARPNESS_COEFFICIENT);
    //const std::vector<TrackPure> restored = restoreTracks(dimensions, grid, responses);
    const std::vector<TrackPure> restored = {
      TrackPure(-0.165431,-0.534174,0.001630,0.025463),
TrackPure(0.724693,-0.455946,0.001833,0.030556),
TrackPure(-0.833663,-2.848085,0.006111,0.020778),
TrackPure(2.499147,3.676024,-0.039111,-0.053778),
TrackPure(6.844856,-7.823867,-0.044000,0.053167),
TrackPure(0.492626,0.851190,-0.026287,-0.036529),
TrackPure(-0.052595,1.639652,-0.029333,0.026889),
TrackPure(-0.600278,0.246349,-0.023222,-0.003056),
TrackPure(4.527494,4.731970,-0.031920,-0.029644),
TrackPure(-0.167654,-0.524886,-0.006722,0.012222),
TrackPure(1.887641,-2.141095,-0.011000,0.018333),
TrackPure(-1.039941,0.494333,-0.023222,0.008556),
TrackPure(0.559952,0.403495,-0.014667,-0.015889),
TrackPure(0.219124,0.075927,-0.019556,-0.005500),
TrackPure(1.769952,3.032667,-0.014667,-0.020778),
TrackPure(-0.391505,-0.427500,-0.015889,0.000000),
TrackPure(0.148496,-0.770665,-0.030556,-0.032389),
TrackPure(-0.212862,0.491043,-0.012833,-0.001222),
TrackPure(0.200427,-0.563865,-0.009778,0.017111),
TrackPure(0.037950,0.957230,-0.022611,0.041556),
TrackPure(4.096792,1.805796,-0.025667,-0.009778),
TrackPure(-1.381505,-0.344802,-0.015889,-0.006111),
TrackPure(-9.061743,9.067835,0.055000,-0.053778),
TrackPure(-3.637729,2.170656,0.025056,-0.011000),
TrackPure(0.017888,0.077604,0.035444,0.040333),
TrackPure(-0.783463,-1.195551,0.009778,0.018333),
TrackPure(-0.592167,0.129755,0.040944,-0.035444),
TrackPure(-0.150963,-0.107989,0.009778,0.004889),
TrackPure(2.159449,-4.172184,0.018333,-0.034833),
TrackPure(-1.183094,4.290385,0.009778,-0.027500),
TrackPure(-2.899563,2.596674,0.034222,-0.028111),
TrackPure(-0.035078,0.388094,0.001222,-0.009778),
TrackPure(-0.252186,-0.679484,-0.004889,-0.008556),
TrackPure(-0.410234,3.380040,0.003667,-0.021389),
TrackPure(-0.630173,4.092962,0.003667,-0.028722),
TrackPure(-4.970255,4.991308,0.021389,0.036667),
TrackPure(-6.879070,-0.610328,0.047056,0.006111),
TrackPure(-2.230091,-1.255726,0.013444,0.007333),
TrackPure(0.871625,-0.825413,0.014667,-0.014056),
TrackPure(0.606048,0.433712,0.014667,0.006722),
TrackPure(-1.765546,2.238376,0.008556,-0.014667),
TrackPure(0.065018,2.736236,-0.001222,-0.020778),
TrackPure(-2.898463,7.849551,0.009778,-0.051944),
TrackPure(-3.383892,3.363814,0.023222,-0.022000),
TrackPure(-2.143952,-6.410260,0.014667,0.040333),
TrackPure(-4.835853,2.548658,0.034833,-0.019556),
TrackPure(0.637172,0.016051,0.015889,-0.004278),
TrackPure(-0.657087,0.168094,0.014056,-0.009778),
TrackPure(-3.253302,-0.547974,0.020778,0.004889),
TrackPure(-0.117578,3.567603,0.001222,-0.022611),
TrackPure(-1.858657,-1.645476,0.019556,-0.012833),
TrackPure(0.259671,-1.051765,0.006111,-0.012222),
TrackPure(-1.504682,1.694652,-0.026889,0.026889),
TrackPure(-1.183212,-4.492629,0.008556,0.028722),
TrackPure(-1.980893,8.154392,-0.019556,-0.061111),
TrackPure(0.165293,0.900532,-0.027539,-0.054053),
TrackPure(1.165318,0.429340,-0.026889,-0.050722),
TrackPure(-0.205825,0.693976,-0.003667,0.012833),
TrackPure(-0.427500,2.975743,0.000000,0.040333),
TrackPure(-0.202339,2.985390,0.011000,0.046444),
TrackPure(0.207500,0.533483,0.000000,-0.047667),
TrackPure(-1.946081,-1.663764,0.018333,-0.020778),
TrackPure(-0.171390,2.123086,-0.011611,0.046444),
TrackPure(-1.731635,0.789782,0.041556,-0.035444),
TrackPure(-2.099526,-0.728468,-0.047322,-0.016238),
TrackPure(-0.644308,-0.301642,-0.049642,-0.015658),
TrackPure(2.020097,1.569173,0.051266,0.032940),
TrackPure(-0.084141,-0.290761,0.016934,-0.015078),
TrackPure(-0.646586,-1.303157,0.018094,-0.012990),
TrackPure(2.368541,-2.016359,-0.034551,0.069679),
TrackPure(7.988331,14.772964,-0.052403,-0.097896),
TrackPure(-0.988105,-10.670357,0.007973,0.071760),
TrackPure(-0.077545,-0.048637,0.009790,0.054131),
TrackPure(11.526244,2.845637,-0.080332,-0.019003),
TrackPure(0.010176,-0.304534,0.008638,0.033400),
TrackPure(0.218903,0.256615,-0.071694,0.032824),
TrackPure(-0.732693,2.385308,-0.022746,0.060177),
TrackPure(8.396873,14.894740,-0.053555,-0.098472),
TrackPure(12.986088,12.774712,-0.088970,-0.085803),
TrackPure(2.854262,6.157126,-0.044840,-0.067948),
TrackPure(10.562625,7.095935,-0.071918,-0.045819),
TrackPure(0.613284,-1.019865,0.033400,-0.059889),
TrackPure(4.716474,14.973269,-0.032248,-0.099624),
TrackPure(-36.836344,-54.221606,0.137237,0.178785),
TrackPure(2.625782,0.265938,0.052115,0.008350),
TrackPure(-13.530423,-4.225531,0.091274,0.028217),
TrackPure(-13.251814,-6.930889,0.088298,0.047028),
TrackPure(-14.653860,-8.713051,0.098184,0.058738),
TrackPure(-1.981991,10.945631,0.014238,-0.074608),
TrackPure(-6.955854,-13.839204,0.046932,0.093865),
TrackPure(0.092619,-0.091360,-0.010722,0.077195),
TrackPure(15.977303,6.527184,-0.106904,-0.044343),
TrackPure(-3.679033,1.406350,-0.072192,0.036811),
TrackPure(1.206804,-0.589380,-0.007386,0.069094),
TrackPure(12.458194,20.190356,-0.083590,-0.135060),
TrackPure(17.131277,18.265054,-0.116150,-0.122225),
TrackPure(-7.990712,-17.288339,0.052893,0.116150),
TrackPure(-0.581879,1.629791,0.032165,-0.075408),
TrackPure(-16.400900,-5.654964,0.109803,0.037820),
TrackPure(7.611842,-5.549537,0.070762,-0.051106),
TrackPure(4.573671,-0.189988,0.071477,-0.001072),
TrackPure(-13.666172,-19.623107,0.092046,0.131315),
TrackPure(20.195309,-0.702565,-0.135849,0.004800),
TrackPure(0.916757,-1.932385,-0.063434,0.107349),
TrackPure(23.100019,16.271225,-0.156158,-0.109385),
TrackPure(20.692359,-21.930301,-0.139559,0.147409),
TrackPure(2.455153,1.779030,-0.096260,-0.065652),
TrackPure(23.922271,19.943872,-0.161574,-0.133756),
TrackPure(7.524057,-3.733391,-0.070975,0.102913),
TrackPure(-3.601670,-20.348372,0.025295,0.132581),
TrackPure(-0.359639,-1.547600,-0.000887,0.062103),
TrackPure(-4.520875,0.279502,-0.028390,-0.100252),
TrackPure(-11.178661,20.403732,0.136479,-0.122136),
TrackPure(-23.990892,0.923478,0.161574,-0.005252),
TrackPure(-4.684758,-22.458408,0.030164,-0.029277),
TrackPure(-23.902246,-26.056578,0.165486,0.173623),
TrackPure(3.037091,0.965018,0.031051,0.111785),
TrackPure(-0.047885,0.016906,0.002218,-0.085613),
TrackPure(-22.307822,-2.841702,0.150291,0.019120),
TrackPure(-0.465529,0.141099,0.034600,-0.089606),
TrackPure(-1.288969,-0.669439,0.091824,0.024841),
TrackPure(-25.613744,22.494300,0.173140,-0.150026),
TrackPure(-0.330674,-0.540071,-0.125033,0.118067),
TrackPure(28.676467,-9.374051,-0.234300,0.030800),
TrackPure(2.018828,-0.368522,-0.118800,0.007700),
TrackPure(37.874507,31.717114,-0.258133,-0.213767),
TrackPure(20.589643,40.899394,-0.137867,-0.274633),
TrackPure(0.506716,0.576872,-0.106333,-0.100467),
TrackPure(37.544645,18.368679,-0.253000,-0.123200),
TrackPure(32.609460,14.656933,-0.218900,-0.099000),
TrackPure(4.440427,-6.610053,-0.067100,0.141900),
TrackPure(30.093388,4.431808,-0.202189,-0.030287),
TrackPure(1.727103,0.969449,0.096800,0.038500),
TrackPure(-1.073678,-0.146504,0.038867,0.129433),
TrackPure(-0.172448,-1.949793,-0.003300,-0.117700),
TrackPure(2.229693,6.910371,0.036667,0.118433),
TrackPure(1.279446,15.301137,-0.017600,-0.177100),
TrackPure(-1.892889,-3.940022,0.077000,0.158400),
TrackPure(-0.074300,0.652029,0.112200,-0.128700),
TrackPure(-0.164858,-3.199776,0.044000,0.141900),
TrackPure(3.675124,35.912666,-0.023917,-0.240932),
TrackPure(-6.754752,-6.958700,0.106839,0.152060),
TrackPure(-3.588334,-6.415552,0.052800,0.140800),
TrackPure(30.528973,-1.634909,-0.275000,0.083600),
TrackPure(-1.870160,-7.437885,-0.025300,-0.122100),
TrackPure(-0.280492,0.402815,0.083600,-0.019800),
TrackPure(-2.948244,0.109569,0.156200,-0.004400),
TrackPure(1.504487,-6.909889,0.033000,-0.111100),
TrackPure(-1.106601,-1.864033,0.048400,0.095700),
TrackPure(-2.861106,0.170776,0.146300,-0.017600),
TrackPure(-0.049587,-1.732195,-0.019800,0.145200),
TrackPure(-0.640831,1.091324,0.161700,0.013200),
TrackPure(-26.838635,-16.446459,0.297367,0.193233),
TrackPure(0.035763,0.252447,0.004400,-0.155100),
TrackPure(12.394530,11.257173,0.102300,0.091300),
TrackPure(4.231342,-25.140105,-0.011000,0.283800),
TrackPure(-4.105742,-7.473882,-0.030800,-0.060500),
TrackPure(2.158608,0.666533,0.035200,0.113300),
TrackPure(-0.326206,-0.289380,0.126368,0.204904),
TrackPure(-1.365892,-4.996275,0.074067,0.273533),
TrackPure(-12.418491,-12.155778,-0.112200,-0.105600),
TrackPure(3.377636,1.018150,-0.196900,-0.059400),
TrackPure(-2.506823,0.277500,-0.179300,0.000000),
TrackPure(0.412667,-0.130116,0.025300,0.165000),
TrackPure(3.584329,3.487912,-0.182600,-0.228800),
TrackPure(0.150382,0.263273,-0.091300,-0.229900),
TrackPure(-16.263542,-2.577674,-0.100833,-0.140067),
TrackPure(2.370146,-3.638643,-0.117333,0.202767),
TrackPure(3.094377,-2.787710,0.198000,-0.145200),
TrackPure(-3.098406,-2.669933,-0.190667,-0.159133),
TrackPure(0.379243,0.490839,-0.171554,-0.275747),
TrackPure(0.198194,-0.450481,0.006600,0.209000),
TrackPure(5.640323,-1.661931,-0.312400,0.099000),
TrackPure(-8.496965,2.264748,0.343567,-0.086533),
TrackPure(-0.308862,0.423657,0.085800,-0.231000),
TrackPure(49.081604,26.661449,-0.329267,-0.178933),
TrackPure(0.229394,0.366971,-0.014300,-0.229900),
TrackPure(38.849230,-23.703384,-0.261076,0.158980),
TrackPure(-21.776714,24.743602,0.146300,-0.166100),
TrackPure(-8.785684,23.176418,0.058667,-0.155833),
TrackPure(-16.713184,32.336581,0.112200,-0.217800),
TrackPure(0.849122,11.044239,0.015400,0.182600),
TrackPure(5.395997,3.159342,-0.290400,-0.182600),
TrackPure(3.297486,-5.186870,-0.189933,0.290767),
TrackPure(4.949744,-1.652394,0.306167,-0.099733),
TrackPure(-5.115995,-27.323658,0.034100,0.183700),
TrackPure(3.045056,0.187178,0.184957,0.009608),
TrackPure(37.197836,-2.420534,-0.248967,0.016133),
TrackPure(19.317051,-30.201767,-0.130533,0.202767),
TrackPure(-15.893759,-0.167385,0.106700,0.001100),
TrackPure(-18.981554,-9.430778,0.127600,0.063800),
TrackPure(-14.298283,19.038544,0.095700,-0.127600),
TrackPure(-4.785262,19.038544,0.031900,-0.127600),
TrackPure(1.332231,15.436795,-0.008800,-0.103400),
TrackPure(5.350798,2.868001,0.325233,0.173067),
TrackPure(15.072437,-6.293180,-0.101200,0.042900),
TrackPure(4.976574,-4.658182,0.289470,-0.281380),
TrackPure(-6.857347,-17.258570,-0.109267,-0.279033),
TrackPure(-24.123863,23.997667,0.278667,-0.276833),
TrackPure(9.570179,21.410618,-0.063800,-0.144100),
TrackPure(19.983390,4.484425,-0.231367,-0.051333),
TrackPure(2.911786,-10.191986,-0.019800,0.068200),
TrackPure(8.487574,20.326331,-0.057200,-0.136400),
TrackPure(16.360153,7.584720,0.266200,0.123200),
TrackPure(8.901758,10.797836,-0.058300,-0.073700),
TrackPure(-3.828990,17.860831,-0.061600,0.292600),
TrackPure(-8.221821,-15.868757,0.053900,0.106700),
TrackPure(17.357197,2.015764,-0.138600,0.004400),
TrackPure(12.755681,5.877019,-0.042900,-0.074800),
TrackPure(10.987835,7.878337,-0.073700,-0.052800),
TrackPure(9.438120,14.416272,-0.109267,-0.167933),
TrackPure(4.495241,-26.350682,0.037400,-0.220000),
TrackPure(6.793233,-3.041998,-0.046200,0.020900),
TrackPure(-7.777438,11.385901,0.055733,-0.078467),
TrackPure(-6.145712,3.045298,0.040700,-0.020900),
TrackPure(5.939160,2.754397,-0.332567,-0.154733),
TrackPure(5.403076,-0.180868,-0.036300,0.002200),
TrackPure(-3.767535,2.923915,0.212890,-0.164959),
TrackPure(0.040069,3.418272,0.002200,-0.136400),
TrackPure(-5.751077,-1.148178,0.084700,0.112200),
TrackPure(-5.867169,-1.860613,0.050600,-0.077000),
TrackPure(0.949327,-3.119012,-0.057200,0.173800),
TrackPure(4.321601,2.451286,-0.048400,-0.028600),
TrackPure(1.494062,-2.742371,-0.085067,0.152533),
TrackPure(2.507197,-0.681496,-0.138600,0.041800),
TrackPure(-4.524486,-1.871109,0.028600,0.013200),
TrackPure(-4.249487,2.560645,0.028600,-0.017600),
TrackPure(0.164058,1.562669,0.009900,-0.115500),
TrackPure(-4.360999,-0.125756,0.030067,0.000733),
TrackPure(2.117699,-0.947521,-0.117700,0.053900),
TrackPure(3.352821,-4.260916,-0.022000,0.028600),
TrackPure(-2.201812,3.012968,0.024839,-0.034695),
TrackPure(-1.390780,-1.107852,0.086738,0.055591),
TrackPure(-0.498696,1.867587,0.021290,-0.074122),
TrackPure(-1.555781,0.542163,0.086738,-0.030753),
TrackPure(-1.287224,-0.301766,0.072019,0.020239),
TrackPure(-0.266675,-0.189680,0.110394,0.064265),
TrackPure(-1.561136,-0.520668,-0.097778,-0.033118),
TrackPure(1.272441,-0.413658,-0.075305,0.025233),
TrackPure(0.046270,-1.509715,0.000789,0.060717),
TrackPure(0.000259,-0.254036,-0.014588,0.087527),
TrackPure(-2.186634,-0.442207,-0.132867,-0.026810),
TrackPure(-2.243460,0.417207,-0.140753,0.026810),
TrackPure(1.391867,0.331194,0.082796,0.021290),
TrackPure(-0.774346,-0.318506,0.042186,0.013405),
TrackPure(-0.165809,0.850889,0.009462,-0.032330),
TrackPure(-3.145746,-0.515150,0.019713,0.003943),
TrackPure(0.841119,-0.264758,0.049677,-0.013011),
TrackPure(3.136550,3.700749,0.071362,0.085950),
TrackPure(-0.219700,-0.292908,0.044946,0.103297),
TrackPure(-0.681549,-1.473896,0.034301,0.083190),
TrackPure(8.926236,0.081992,0.145244,0.000530),
TrackPure(6.560190,5.521943,0.105137,0.089630),
TrackPure(0.254969,-2.892022,-0.001577,0.018925),
TrackPure(-4.312602,-0.968746,-0.098172,-0.020502),
TrackPure(0.924542,-1.938962,-0.019618,0.111300),
TrackPure(-1.390828,-6.714734,-0.026824,-0.158543),
TrackPure(2.963898,27.416034,0.024822,0.228606),
TrackPure(-4.174509,3.442500,-0.095938,0.082270),
TrackPure(-3.113779,5.530390,-0.050466,0.089498),
TrackPure(0.968807,-0.392383,-0.038638,0.016165),
TrackPure(0.018426,0.052016,-0.057563,0.017348),
TrackPure(1.781609,0.407093,0.035200,0.009900),
TrackPure(0.483773,0.128895,0.029700,0.006600),
TrackPure(0.071969,-2.010801,0.002200,-0.048400),
TrackPure(0.012501,2.316494,0.000000,-0.016500),
TrackPure(-0.200525,-0.190243,0.033000,0.046200),
TrackPure(1.210950,17.381860,0.008800,0.144100),
TrackPure(-23.121824,15.182425,-0.194700,0.126500),
TrackPure(-0.747726,3.203334,-0.013200,0.052800),
TrackPure(-14.662915,4.906295,-0.122496,0.042207),
TrackPure(-0.280758,-0.029043,0.004400,0.034100),
TrackPure(4.522327,4.593067,-0.029700,-0.030800),
TrackPure(0.422476,-0.089555,-0.024933,0.001833),
TrackPure(-7.829629,9.968545,-0.064533,0.082867),
TrackPure(0.055121,0.036310,-0.026400,-0.013200),
TrackPure(-0.250475,0.075160,-0.033000,0.011000),
TrackPure(0.686586,-0.270138,-0.028600,0.007700),
TrackPure(4.169786,-24.900918,0.034469,-0.208221),
TrackPure(3.640465,-4.043312,-0.205333,0.229167),
TrackPure(2.300102,4.741176,-0.129067,-0.266933),
TrackPure(-2.828432,-10.997256,-0.047300,-0.172700),
TrackPure(4.016945,27.774064,0.030800,0.229900),
TrackPure(-2.703068,-0.387945,-0.140800,0.024200),
TrackPure(-1.018063,-3.144673,0.040700,0.005500),
TrackPure(0.072422,-0.636116,-0.063800,0.215600),
TrackPure(0.669227,-0.801817,-0.280867,0.275733),
TrackPure(-2.980438,-1.775827,-0.182600,-0.105600),
TrackPure(2.069399,7.077019,0.048400,0.167200),
TrackPure(-1.132227,1.879552,-0.038638,0.056380),
TrackPure(-8.025700,17.458612,-0.093151,0.212325),
TrackPure(0.468573,4.077891,0.029700,0.249700),
TrackPure(10.550535,8.785705,-0.035484,-0.037455),
TrackPure(-2.967257,0.562038,-0.143118,0.028781),
TrackPure(-4.758210,7.963671,-0.110657,0.186225),
TrackPure(11.530904,7.462550,-0.031655,-0.044701),
TrackPure(-0.386274,-0.218344,-0.007333,-0.004889),
TrackPure(0.969274,-3.569446,0.007333,0.038500),
TrackPure(-0.297345,-0.592063,-0.002444,-0.007333),
TrackPure(0.191837,3.381765,0.018700,0.016500),
TrackPure(3.939093,6.648721,-0.143000,-0.257400),
TrackPure(-1.482916,1.050164,0.010389,-0.007944),
TrackPure(5.841376,-7.329509,-0.001222,-0.023833),
TrackPure(-55.336305,20.545673,0.371518,-0.138228),
TrackPure(-52.242022,-9.833162,0.350319,0.066860),
TrackPure(32.248200,41.558070,-0.219267,-0.281233),
TrackPure(-30.086176,-33.870774,0.202400,0.229900),
TrackPure(-3.711500,-15.760549,0.019775,0.113052),
TrackPure(-17.022546,-15.700141,0.108716,0.101962),
TrackPure(-2.544507,-17.308820,0.015400,0.118800),
TrackPure(3.968344,-24.817131,0.279767,0.011733),
TrackPure(-10.250158,-27.689461,0.069667,0.187000),
TrackPure(-10.190199,-27.767333,0.069300,0.191400),
TrackPure(-15.121882,-14.719613,0.101933,0.098267),
TrackPure(47.465519,30.896901,-0.319367,-0.207533),
TrackPure(48.514916,-0.052502,-0.325600,0.000000),
TrackPure(3.532186,9.206820,0.019800,-0.220000),
TrackPure(7.255778,-2.126105,-0.294433,0.087267),
TrackPure(-4.442549,22.424918,0.031167,-0.151067),
TrackPure(-4.285082,-2.159868,0.236500,0.119900),
TrackPure(0.848286,-0.042813,-0.312400,0.019800) 
    };
    
    /*{
      std::ofstream myfile;
      myfile.open ("tracks.csv");
      myfile << "x0,y0,dx,dy,reps" << std::endl;
      for (const TrackPure& track: restored)
      {
        myfile << track.xOnZ0 << "," << track.yOnZ0 << "," 
          << track.dxOverDz << "," << track.dyOverDz << ","
          << calculateResponses({track}, hits, RETINA_SHARPNESS_COEFFICIENT)[0]<< std::endl;
      }
      myfile.close();
    }*/
    std::vector<TrackPure> maxims;

    {
      std::ofstream myfile;
      myfile.open ("maxims.csv");
      myfile << "x0,y0,dx,dy,reps" << std::endl;
      std::default_random_engine generator;
      std::normal_distribution<double> dx(0, 0.01);
      std::normal_distribution<double> x0(0, 0.3);
      std::uniform_int_distribution<int> distribution(0, hits.size());  

      int cnt = 0, nw = 0;
      while (cnt < 10000 && maxims.size() < 2000)
      //for (size_t i = 0 ; i < )
      {
        TrackPure track;
        {
          /*int a, b;
          do
          {
            a = distribution(generator);
            b = distribution(generator);            
          } while (hits[a].sensorId == hits[b].sensorId);
          track = TrackPure(hits[a], hits[b]);*/
          track = makeLocalMaximum(
            TrackPure(
              x0(generator),
              x0(generator),
              dx(generator),
              dx(generator)
            ),
            hits, 
            RETINA_SHARPNESS_COEFFICIENT
            );
        }
        if ((calculateResponses({track}, hits, RETINA_SHARPNESS_COEFFICIENT)[0] > 2.1) && 
        std::all_of(maxims.begin(), maxims.end(), 
          [&](TrackPure& old)
          {
            if (fabs(old.xOnZ0 - track.xOnZ0) > 0.1)
              return true;
            if (fabs(old.yOnZ0 - track.yOnZ0) > 0.1)
              return true;
            if (fabs(old.dxOverDz - track.dxOverDz) > 0.01)
              return true;
            if (fabs(old.dyOverDz - track.dyOverDz) > 0.01)
              return true;
            return false;
          }) 
          )
        {
          maxims.push_back(track);
          nw = 0;
          myfile << track.xOnZ0 << "," << track.yOnZ0 << "," 
            << track.dxOverDz << "," << track.dyOverDz << ","
            << calculateResponses({track}, hits, RETINA_SHARPNESS_COEFFICIENT)[0]<< std::endl;
        }
        else
        {
          nw++;
        }
        if (++cnt % 10 == 0)
        {
          std::cerr << cnt << " tracks procceed:" << maxims.size() <<  std::endl;
        }
 
      }
    
      myfile.close();
    }
    
    /*int cnt = 0;
    for (const TrackPure& track: restored)
    {
      auto nw = makeLocalMaximum(track, hits, RETINA_SHARPNESS_COEFFICIENT);
      if (calculateResponses({nw}, hits, RETINA_SHARPNESS_COEFFICIENT)[0] > 2)
        maxims.push_back(nw);
      if (++cnt % 100 == 0)
      {
        std::cerr << cnt << " tracks procceed" << std::endl;
      }
    }*/
    const std::vector<Track> tracksWithHits = findHits(maxims, hits);
    auto answer = putTracksInOutputFormat(hits, tracksWithHits);
    output[i] = answer;
    //printSolution(tracksWithHits, hits, DEBUG);
  }
  return 0;
}
