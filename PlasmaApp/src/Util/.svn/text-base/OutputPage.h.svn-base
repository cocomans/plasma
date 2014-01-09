#ifndef OUTPUT_PAGE_H
#define OUTPUT_PAGE_H
#include <string.h>
#include <sys/stat.h>    /* Fix up for Windows - inc mode_t */
#include <iostream>

class OutputPage
{
public:

	OutputPage(){
	iline = -1;
	nlines = 0;
	nlines_allocated = 128;
	lines = (std::string**)malloc(nlines_allocated*sizeof(std::string*));}

	std::string &nextline()
	{

		iline += 1;

		return getline(iline);

	}
	std::string &prevline()
	{
		iline -= 1;
		return getline(iline);
	}

	std::string &getline(int i)
	{
		if(i >= nlines)
		{
			lines[i] = new std::string("");
			nlines = i+1;
		}

		if(nlines >= nlines_allocated)
		{
			// Need to allocate more lines
			nlines_allocated = (nlines+128)/128;
			nlines_allocated *= 128;
			std::string** lines2 = (std::string**)malloc(nlines_allocated*sizeof(std::string*));
			memcpy(lines2,lines,nlines*sizeof(std::string*));
			free(lines);
			lines = lines2;
		}

		i = std::max(i,0);

		return *(lines[i]);

	}

	void writepage(const char* filename)
	{
		FILE* fp;
		iline = -1;

		fp = fopen(filename,"w");

		for(int i=0;i<nlines;i++)
		{
			fprintf(fp,"%s\n",nextline().c_str());
		}

		fclose(fp);
	}

protected:
	std::string** lines;
	int nlines;
	int nlines_allocated;
	int iline;

};

#endif /* OUTPUT_PAGE_H */
