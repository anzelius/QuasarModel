import os


# lazy quick fix so program finds built in fits files when run from both cl or main
two_sources = 'fits/two_sources/J0000+4054_X_2014_06_09_pet_map.fits' \
    if os.path.exists('fits/two_sources/J0000+4054_X_2014_06_09_pet_map.fits') \
    else 'quasarModel/fits/two_sources/J0000+4054_X_2014_06_09_pet_map.fits'
one_source1 = 'fits/smeared_source/J0029-0113_S_2004_04_30_yyk_map.fits' \
    if os.path.exists('fits/smeared_source/J0029-0113_S_2004_04_30_yyk_map.fits') \
    else 'quasarModel/fits/smeared_source/J0029-0113_S_2004_04_30_yyk_map.fits'
one_source2 = 'fits/small_oval_no_tilt/J0240+0731_C_2016_07_06_pet_map.fits' \
    if os.path.exists('fits/small_oval_no_tilt/J0240+0731_C_2016_07_06_pet_map.fits') \
    else 'quasarModel/fits/small_oval_no_tilt/J0240+0731_C_2016_07_06_pet_map.fits'
one_source3 = 'fits/thin_oval_tilt/J0226-4255_C_2014_02_23_pet_map.fits' \
    if os.path.exists('fits/thin_oval_tilt/J0226-4255_C_2014_02_23_pet_map.fits') \
    else 'quasarModel/fits/thin_oval_tilt/J0226-4255_C_2014_02_23_pet_map.fits'
one_source4 = 'fits/oval_tilt/J0211-0145_C_2019_07_24_pet_map.fits' \
    if os.path.exists('fits/oval_tilt/J0211-0145_C_2019_07_24_pet_map.fits') \
    else 'quasarModel/fits/oval_tilt/J0211-0145_C_2019_07_24_pet_map.fits'
one_source5 = 'fits/small_oval_tilt/J0233+3115_C_2019_08_01_pet_map.fits' \
    if os.path.exists('fits/small_oval_tilt/J0233+3115_C_2019_08_01_pet_map.fits') \
    else 'quasarModel/fits/small_oval_tilt/J0233+3115_C_2019_08_01_pet_map.fits'
one_source6 = 'fits/oval_tilt_noisy/J0236-1445_C_2015_08_31_pet_map.fits' \
    if os.path.exists('fits/oval_tilt_noisy/J0236-1445_C_2015_08_31_pet_map.fits') \
    else 'quasarModel/fits/oval_tilt_noisy/J0236-1445_C_2015_08_31_pet_map.fits'

help_analytical = 'Run program in terminal with analytical source inputs. Input parameters for analytical source will be ' \
          'prompted in terminal.'
help_real = 'Run program in terminal with real sources in fits format. Provide the filepath for the desired fits file' \
          'or use one of the seven test files provided in the program by using corresponding number as input. ' \
          'Description of test files: \n' \
          '1: One oval source which is a bit smeared. \n ' \
          '2: One small oval source with no tilt. \n ' \
          '3: One thin oval tilted source \n' \
          '4: One a bit thicker oval tilted source \n' \
          '5: One small oval tilted source \n' \
          '6: One tilted oval and a bit noisy source \n' \
          '7: Two sources'
help_directory = 'Run program on several sources by providing a directory containing source directories with fits files'

stars = '*********************************************'
line = '-----------------------------------------------------------------------------------------------'

icon_path = 'images/lovecat.ico' if os.path.exists('images/lovecat.ico') else 'quasarModel/images/lovecat.ico'