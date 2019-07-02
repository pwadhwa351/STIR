//
//
/*
    Copyright (C) 2000 PARAPET partners
    Copyright (C) 2000- 2008, Hammersmith Imanet Ltd
    This file is part of STIR.

    This file is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.

    This file is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

	See STIR/LICENSE.txt for details
*/

#ifndef __stir_config_gcc_H__
#define __stir_config_gcc_H__

/*!
  \file 
  \ingroup buildblock 
  \brief configuration for gcc

  \author Kris Thielemans
  \author PARAPET project




 This include file defines a few macros and en/disables pragmas
 specific to gcc.

 It is included by sitr/common.h. You should never include it directly.
*/

#if defined __GNUC__
# if __GNUC__ == 2 && __GNUC_MINOR__ <= 8
#  define STIR_NO_NAMESPACES
#  define STIR_NO_AUTO_PTR
# endif
#endif

#endif 
