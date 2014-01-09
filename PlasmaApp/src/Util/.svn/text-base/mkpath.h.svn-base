/*
@(#)File:           $RCSfile: mkpath.h,v $
@(#)Version:        $Revision: 1.13 $
@(#)Last changed:   $Date: 2012/07/15 00:40:37 $
@(#)Purpose:        Create all directories in path
@(#)Author:         J Leffler
@(#)Copyright:      (C) JLSS 1990-91,1997-98,2001,2005,2008,2012
*/

#ifndef MAKE_PATH_H
#define MAKE_PATH_H
#include <errno.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif /* HAVE_UNISTD_H */
#include <string.h>
#include <sys/stat.h>    /* Fix up for Windows - inc mode_t */

typedef struct stat Stat;

static int do_mkdir(const char *path, mode_t mode);
/**
** mkpath - ensure all directories in path exist
** Algorithm takes the pessimistic view and works top-down to ensure
** each directory in path exists, rather than optimistically creating
** the last element and working backwards.
*/
int mkpath(const char *path, mode_t mode);
#endif /* MAKE_PATH_H */
