//
// $Id$
//
/*!

  \file
  \ingroup OSMAPOSL

  \brief This is a program to compute the 'sensitivity' image (detection probabilities per voxel). 
  When no input attenuation file is specified, attenuation factors are supposed to be 1, and similar for the normalisation factors.


  \author Matthew Jacobson
  \author Kris Thielemans
  \author Alexey Zverovich
  \author PARAPET project

  $Date$
  $Revision$
*/

/* 

   if TEST is #defined, results of different ring differences (i.e. the profiles   through the centres of all planes) are
   stored as separate files (seg_x.dat) and the sensitivity image is
   not computed.
   if TEST is not defined, output is written to a user-specified output file

 
   TODO:
   - implement the Normalisation class
   */


#include "ProjData.h"
#include "ProjDataFromStream.h"
#include "RelatedViewgrams.h"
#include "ProjDataInfoCylindrical.h"
#include "VoxelsOnCartesianGrid.h"
#include "utilities.h"
#include "interfile.h"

#include "ArrayFunction.h"
#include "recon_array_functions.h"

#include "LogLikBased/common.h"

#include "recon_buildblock/distributable.h"
#if 0
#ifndef USE_PMRT
#include "recon_buildblock/ForwardProjectorByBinUsingRayTracing.h"
#include "recon_buildblock/BackProjectorByBinUsingInterpolation.h"
#else
#include "recon_buildblock/ForwardProjectorByBinUsingProjMatrixByBin.h"
#include "recon_buildblock/BackProjectorByBinUsingProjMatrixByBin.h"
#include "recon_buildblock/ProjMatrixByBinUsingRayTracing.h"
#endif
#else
#include "recon_buildblock/ForwardProjectorByBin.h"
#include "recon_buildblock/BackProjectorByBin.h"
#endif

#include "OSMAPOSL/OSMAPOSLParameters.h"

#ifndef TOMO_NO_NAMESPACES
using std::endl;
#endif


START_NAMESPACE_TOMO

class sensparameters{
  public:
 int limit_segments;
  bool zero_seg0_end_planes;
  int rim_truncation_sino;
} globals;



// A do-nothing class for normalisation
class Normalisation
{
public:
  virtual void apply(RelatedViewgrams<float>&) const {}
};

// KT 01/12/98 function which constructs ProjData (nothing else for the moment)
ProjData * ask_parameters();



void
 compute_sensitivity_image(DiscretisedDensity<3,float>& result,
                           shared_ptr<ProjData> const& proj_data_ptr,
			   shared_ptr<DiscretisedDensity<3,float> > const& attenuation_image_ptr,
			   const bool do_attenuation,
                           shared_ptr<ProjData> const& atten_proj_data_ptr, 
			   shared_ptr<ProjData> const& norm_proj_data_ptr, 
			   //const Normalisation& normalisation,
			   const sensparameters &globals)
{





#ifdef TEST

  for (int segment_num = 0; segment_num <= globals.limit_segments; segment_num++)
  {
    int min_segment = segment_num;
    int max_segment = segment_num;
    
    result.fill(0);
      
#else

    int min_segment = 0;
    int max_segment = globals.limit_segments;

#endif

    //warning("Currently ignoring attenuation factors (but not the attenuation image)\n");
    // TODO norm
    distributable_compute_sensitivity_image(result,
					    proj_data_ptr->get_proj_data_info_ptr()->clone(),
                                            do_attenuation ? attenuation_image_ptr.get() : NULL,
					    0,
					    1,
					    min_segment,
					    max_segment,
					    globals.zero_seg0_end_planes,
					    norm_proj_data_ptr);

#ifdef TEST

    char fname[20];
    sprintf(fname, "seg_%d.prof", segment_num);
    cerr << "Writing horizontal profiles to " << fname << endl;
    ofstream profile(fname);
    if (!profile)
    {
      cerr << "Couldn't open " << fname << endl;
    }
    else
    {
      for (int z = result.get_min_z(); z <= result.get_max_z(); z++) 
      { 
	for (int x = result.get_min_x(); x <= result.get_max_x(); x++)
	  profile<<(*result)[z][0][x]<<" ";
	profile << "\n";
      };
    };
  };

#endif

}


ProjData * 
ask_parameters()
{
  
  shared_ptr<ProjDataInfo> proj_data_info_ptr = ProjDataInfo::ask_parameters();

  iostream *out = 0;
  
  return new ProjDataFromStream(
    proj_data_info_ptr, 
    out, 0UL);
}




END_NAMESPACE_TOMO

USING_NAMESPACE_TOMO

#ifdef PARALLEL
int master_main(int argc, char **argv)
#else
int main(int argc, char **argv)
#endif
{
  shared_ptr<ProjData> proj_data_ptr = 0;


#if 0  
  if(argc==2)
  {
    proj_data_ptr = ProjData::read_from_file(argv[1]);
  }
  else
  {
    proj_data_ptr = ask_parameters();
  }
  if(argc!=2) 
  {
    cerr<<"Usage: sensitivity [PSOV-file]\n"
        <<"The PSOV-file will be used to get the scanner, mashing etc. details" 
	<< endl; 
  }
   

  globals.limit_segments=ask_num("Maximum absolute segment number to process: ", 
    0, proj_data_ptr->get_max_segment_num(), proj_data_ptr->get_max_segment_num() );
  
  {
    const ProjDataInfoCylindrical * pdi_cyl_ptr =
      dynamic_cast<const ProjDataInfoCylindrical *>
      ( proj_data_ptr->get_proj_data_info_ptr());
    
    globals.zero_seg0_end_planes =
      pdi_cyl_ptr != NULL &&
      ask("Zero end planes of segment 0 ?", 
           pdi_cyl_ptr->get_min_ring_difference(0) == pdi_cyl_ptr->get_max_ring_difference(0));
  }

  globals.rim_truncation_sino = rim_truncation_sino;
    // ask_num("Number of bins to set to zero ?",0, proj_data_ptr->get_max_tangential_pos_num(), 4);


#ifndef TEST
  char out_filename[max_filename_length];
  ask_filename_with_extension(out_filename,
			      "Output to which file (without extension)?",
			      "");
#endif

#else
  if(argc!=2) 
  {
    cerr<<"Usage: sensitivity OSMAPOSL_par_file\n"
        <<"The par-file will be used to get the scanner, mashing etc. details" 
	<< endl; 
  }
  if (argc>2)
    exit(1);

  OSMAPOSLParameters parameters(argv[1]);
  if (parameters.max_segment_num_to_process==-1)
    parameters.max_segment_num_to_process =
      parameters.proj_data_ptr->get_max_segment_num();

  globals.rim_truncation_sino = rim_truncation_sino;
  string out_filename = parameters.sensitivity_image_filename;
  globals.zero_seg0_end_planes = parameters.zero_seg0_end_planes;
  proj_data_ptr = parameters.proj_data_ptr;
  globals.limit_segments = parameters.max_segment_num_to_process;
  
#endif

  bool do_attenuation;
  shared_ptr<DiscretisedDensity<3,float> >
    attenuation_image_ptr(
    new VoxelsOnCartesianGrid<float>(*proj_data_ptr->get_proj_data_info_ptr()));

  {
    char atten_name[max_filename_length];
    ask_filename_with_extension(atten_name, 
				"Get attenuation image from which file (0 = 0's): ",
				"");    
    
    // KT 18/08/2000 compare whole string instead of only first character
    if(strcmp(atten_name,"0")==0)
    {
      do_attenuation = false;
    }
    else
    {
      do_attenuation = true;

     // Read from file
      // TODO remove unsatisfactory manipulations for handling 'odd-sized' images.
      // TODO add checks on sizes etc.

      shared_ptr<DiscretisedDensity<3,float> > attenuation_density_from_file_ptr =
        DiscretisedDensity<3,float>::read_from_file(atten_name);

      const VoxelsOnCartesianGrid<float> * attenuation_image_from_file_ptr =
         dynamic_cast<const VoxelsOnCartesianGrid<float> *>(attenuation_density_from_file_ptr.get());
      
      if (attenuation_image_from_file_ptr == NULL)
        error("Can only handle VoxelsOnCartesianGrid for the attenuation image\n");

      // Now add to the attenuation_image (which is 0 at this point)
      // This in principle should allow us to read an 'even-sized' image
      // as += just adds the appropriate ranges
      *attenuation_image_ptr += *attenuation_image_from_file_ptr;

      cerr << "WARNING: attenuation image data are supposed to be in units cm^-1\n"
	"Reference: water has mu .096 cm^-1" << endl;
      cerr<< "Max in attenuation image:" 
	  << attenuation_image_ptr->find_max() << endl;
#ifndef NORESCALE
      /*
	cerr << "WARNING: multiplying attenuation image by x-voxel size "
	<< " to correct for scale factor in forward projectors...\n";
      */
      // projectors work in pixel units, so convert attenuation data 
      // from cm^-1 to pixel_units^-1
      const float rescale = 
	dynamic_cast<VoxelsOnCartesianGrid<float> *>(attenuation_image_ptr.get())->
	get_voxel_size().x()/10;
#else
      const float rescale = 
	dynamic_cast<VoxelsOnCartesianGrid<float> *>(attenuation_image_ptr.get())->
	10.F;
#endif
      *attenuation_image_ptr *= rescale;      
    }
  }

  shared_ptr<ProjData> atten_proj_data_ptr = 0;
#if 0
  // TODO
  if (!do_attenuation)
  {
    char atten_name[max_filename_length];
    ask_filename_with_extension(atten_name, 
				"Get attenuation factors from which file (1 = 1's):",
				"");    
    
    if(strcmp(atten_name,"1")!=0)
    {
      atten_proj_data_ptr = ProjData::read_from_file(atten_name);
    }
  }
#endif

  shared_ptr<ProjData> norm_proj_data_ptr = 0;
#if 1
  {
    char norm_name[max_filename_length];
    ask_filename_with_extension(norm_name, 
				"Get normalisation factors from which file (1 = 1's):",
				"");    
    
    if(strcmp(norm_name,"1")!=0)
    {
      norm_proj_data_ptr = ProjData::read_from_file(norm_name);
    }
  }
#endif
#if 0
  // TODO replace
#ifndef USE_PMRT
  shared_ptr<ForwardProjectorByBin> forward_projector_ptr =
    new ForwardProjectorByBinUsingRayTracing(proj_data_ptr->get_proj_data_info_ptr()->clone(), 
                                             attenuation_image_ptr);
  shared_ptr<BackProjectorByBin> back_projector_ptr =
    new BackProjectorByBinUsingInterpolation(proj_data_ptr->get_proj_data_info_ptr()->clone(), 
                                             attenuation_image_ptr);
#else
  shared_ptr<ProjMatrixByBin> PM = 
    new  ProjMatrixByBinUsingRayTracing( attenuation_image_ptr , proj_data_ptr->get_proj_data_info_ptr()->clone()); 	
  ForwardProjectorByBin* forward_projector_ptr =
    new ForwardProjectorByBinUsingProjMatrixByBin(PM); 
  BackProjectorByBin* back_projector_ptr =
    new BackProjectorByBinUsingProjMatrixByBin(PM); 
#endif

#else
  shared_ptr<ForwardProjectorByBin> forward_projector_ptr =
    parameters.projector_pair_ptr->get_forward_projector_sptr();
  shared_ptr<BackProjectorByBin> back_projector_ptr =
    parameters.projector_pair_ptr->get_back_projector_sptr();
  parameters.projector_pair_ptr->set_up(proj_data_ptr->get_proj_data_info_ptr()->clone(), 
					attenuation_image_ptr);
#endif
  set_projectors_and_symmetries(forward_projector_ptr, 
                                back_projector_ptr, 
                                back_projector_ptr->get_symmetries_used()->clone());


  // Initialise the sensitivity image  
  shared_ptr<DiscretisedDensity<3,float> > result_ptr;
#if 0
  result_ptr =
    attenuation_image_ptr->get_empty_discretised_density();
#else
    // TODO somehow get rid of VoxelsOnCartesianGrid
  if(parameters.initial_image_filename=="1")
  {
    result_ptr =
      new VoxelsOnCartesianGrid<float> (*parameters.proj_data_ptr->get_proj_data_info_ptr(),
					parameters.zoom,
					CartesianCoordinate3D<float>(parameters.Zoffset,
								     parameters.Yoffset,
								     parameters.Xoffset),
					parameters.output_image_size);
  }
  else
    {
      result_ptr = 
        DiscretisedDensity<3,float>::read_from_file(parameters.initial_image_filename);
      result_ptr->fill(0);
    }
#endif

  // Compute the sensitivity image  
    compute_sensitivity_image(*result_ptr,
			      proj_data_ptr, 
                              attenuation_image_ptr,  do_attenuation, 
                              atten_proj_data_ptr,
                              norm_proj_data_ptr,
                              //Normalisation (), 
                              globals);
#ifndef TEST
  write_basic_interfile(out_filename, *result_ptr, NumericType::FLOAT);

  cerr << "min and max in image " << result_ptr->find_min() 
       << " " << result_ptr->find_max() << endl;

#endif

  return EXIT_SUCCESS;

}
