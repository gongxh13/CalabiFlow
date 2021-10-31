/*!
*      \file RgbImage.h
*      \brief BitMap Image
*      \date Documented on 10/08/2010
*
*/
#ifndef RGBIMAGE_H
#define RGBIMAGE_H

#include <stdio.h>
#include <assert.h>

#ifndef RGBIMAGE_DONT_USE_OPENGL
#include <windows.h>
#include "GL/gl.h"
#endif

// Include the next line to turn off the routines that use OpenGL
// #define RGBIMAGE_DONT_USE_OPENGL
/*!
 *	\brief RgbImage class
 *
 *  24 bit bmp image class.
 */
class RgbImage
{
public:
	/*!	RgbImage constructor
	*/
	RgbImage();
	/*!	RgbImage constructor
	* \param filename the input file nmae
	*/
	RgbImage( const char* filename );
	/*! RgbImage constructor
	*	\param numRows number of rows
	*   \param numCols number of cols
	*/
	RgbImage( int numRows, int numCols );	// Initialize a blank bitmap of this size.
	/*!	RgbImage destructor
	 *
	 */
	~RgbImage();
	
	/*!
	 *	LoadBmpFile 
	 * \param filename input file name
	 * \return whether the loading is successful
	 */
	bool LoadBmpFile( const char *filename );		// Loads the bitmap from the specified file
	/*!
	 *	WriteBmpFile 
	 * \param filename output file name
	 * \return whether the writing is successful
	 */
	bool WriteBmpFile( const char* filename );		// Write the bitmap to the specified file
#ifndef RGBIMAGE_DONT_USE_OPENGL
	/*! Load the Bmp Image from OpenGL buffer
	*  \return whether the loading is successful
	*/
	bool LoadFromOpenglBuffer();					// Load the bitmap from the current OpenGL buffer
#endif

	/*!	number of rows
	*/
	long GetNumRows() const { return NumRows; }
	/*!	number of columns
	*/
	long GetNumCols() const { return NumCols; }
	// Rows are word aligned
	/*! Get number of bytes per row
	 */	
	long GetNumBytesPerRow() const { return ((3*NumCols+3)>>2)<<2; }	
	/*!	Get image data buffer
	*/
	const void* ImageData() const { return (void*)ImagePtr; }
	/*! Get the pixel at (row, col)
	 * \param row row position
	 * \param col column position
	 * \return array of unsigned char for the pixel color
	 */
	const unsigned char* GetRgbPixel( long row, long col ) const;
	/*! Get the pixel at (row, col)
	 * \param row row position
	 * \param col column position
	 * \return array of unsigned char for the pixel color
	 */
	unsigned char* GetRgbPixel( long row, long col );
	/*! Get the pixel at (row, col)
	 * \param row row position
	 * \param col column position
	 * \param red pointer to the red color
	 * \param green pointer to the green color
	 * \param blue  pointer to the blue  color
	 */
	void GetRgbPixel( long row, long col, float* red, float* green, float* blue ) const;
	/*! Get the pixel at (row, col)
	 * \param row row position
	 * \param col column position
	 * \param red pointer to the red color
	 * \param green pointer to the green color
	 * \param blue  pointer to the blue  color
	 */
	void GetRgbPixel( long row, long col, double* red, double* green, double* blue ) const;

	/*! set the pixel at (row, col)
	 * \param row row position
	 * \param col column position
	 * \param red   the red color floating point
	 * \param green the green color floating point
	 * \param blue  the blue  color floating point
	 */
	void SetRgbPixelf( long row, long col, double red, double green, double blue );

	/*! set the pixel at (row, col)
	 * \param row row position
	 * \param col column position
	 * \param red   the red color unsigned char
	 * \param green the green color unsigned char
	 * \param blue  the blue  color unsigned char
	 */
	void SetRgbPixelc( long row, long col, 
					   unsigned char red, unsigned char green, unsigned char blue );
	/*! get the error code
	 */

	// Error reporting. (errors also print message to stderr)
	int GetErrorCode() const { return ErrorCode; }
	enum {
		NoError = 0,
		OpenError = 1,			// Unable to open file for reading
		FileFormatError = 2,	// Not recognized as a 24 bit BMP file
		MemoryError = 3,		// Unable to allocate memory for image data
		ReadError = 4,			// End of file reached prematurely
		WriteError = 5			// Unable to write out data (or no date to write out)
	};
	bool ImageLoaded() const { return (ImagePtr!=0); }  // Is an image loaded?

	void Reset();			// Frees image data memory

private:
	/*!	
	pixel buffer
	*/
	unsigned char* ImagePtr;	// array of pixel values (integers range 0 to 255)
	/*! number of rows
	*/
	long NumRows;				// number of rows in image
	/*! number of columns
	*/
	long NumCols;		// number of columns in image
	/*! current error code
	*/
	int ErrorCode;				// error code
	
	/*!	read a short integer form a file
	\param infile input file
	\return short integer just read
	*/
	static short readShort( FILE* infile );
	/*!	read a long integer form a file
	\param infile input file
	\return long integer just read
	*/
	static long readLong( FILE* infile );
	/*!	skip several chars in the file
	\param infile input file name
	\param numChars the number of chars to be skipped
	*/
	
	static void skipChars( FILE* infile, int numChars );
	/*!	write a long integer to the file
	\param data the long integter to write
	\param outfile output file
	*/
	static void RgbImage::writeLong( long data, FILE* outfile );
	/*!	write a short integer to the file
	\param data the short integter to write
	\param outfile output file
	*/
	static void RgbImage::writeShort( short data, FILE* outfile );
	
	/*!	convert a double to unsigned char
	* \param x input double number
	* \return output unsigned char
	*/
	static unsigned char doubleToUnsignedChar( double x );

};

/*!	RgbImage constructor
	*/
inline RgbImage::RgbImage()
{ 
	NumRows = 0;
	NumCols = 0;
	ImagePtr = 0;
	ErrorCode = 0;
}
/*!	RgbImage constructor
	* \param filename the input file nmae
	*/
inline RgbImage::RgbImage( const char* filename )
{
	NumRows = 0;
	NumCols = 0;
	ImagePtr = 0;
	ErrorCode = 0;
	LoadBmpFile( filename );
}
/*!	RgbImage destructor
	 *
	 */
inline RgbImage::~RgbImage()
{ 
	delete[] ImagePtr;
}

// Returned value points to three "unsigned char" values for R,G,B
/*! Get the pixel at (row, col)
	 * \param row row position
	 * \param col column position
	 * \return array of unsigned char for the pixel color
	 */
inline const unsigned char* RgbImage::GetRgbPixel( long row, long col ) const
{
	assert ( row<NumRows && col<NumCols );
	const unsigned char* ret = ImagePtr;
	long i = row*GetNumBytesPerRow() + 3*col;
	ret += i;
	return ret;
}
/*! Get the pixel at (row, col)
	 * \param row row position
	 * \param col column position
	 * \return array of unsigned char for the pixel color
	 */
inline unsigned char* RgbImage::GetRgbPixel( long row, long col ) 
{
	assert ( row<NumRows && col<NumCols );
	unsigned char* ret = ImagePtr;
	long i = row*GetNumBytesPerRow() + 3*col;
	ret += i;
	return ret;
}
/*! Get the pixel at (row, col)
	 * \param row row position
	 * \param col column position
	 * \param red pointer to the red color
	 * \param green pointer to the green color
	 * \param blue  pointer to the blue  color
	 */
inline void RgbImage::GetRgbPixel( long row, long col, float* red, float* green, float* blue ) const
{
	assert ( row<NumRows && col<NumCols );
	const unsigned char* thePixel = GetRgbPixel( row, col );
	const float f = 1.0f/255.0f;
	*red = f*(float)(*(thePixel++));
	*green = f*(float)(*(thePixel++));
	*blue = f*(float)(*thePixel);
}
/*! Get the pixel at (row, col)
	 * \param row row position
	 * \param col column position
	 * \param red pointer to the red color
	 * \param green pointer to the green color
	 * \param blue  pointer to the blue  color
	 */
inline void RgbImage::GetRgbPixel( long row, long col, double* red, double* green, double* blue ) const
{
	assert ( row<NumRows && col<NumCols );
	const unsigned char* thePixel = GetRgbPixel( row, col );
	const double f = 1.0/255.0;
	*red = f*(double)(*(thePixel++));
	*green = f*(double)(*(thePixel++));
	*blue = f*(double)(*thePixel);
}
/*! reset the RgbImage object
*/
inline void RgbImage::Reset()
{
	NumRows = 0;
	NumCols = 0;
	delete[] ImagePtr;
	ImagePtr = 0;
	ErrorCode = 0;
}


/*!
*      \file RgbImage.cpp
*      \brief Implementation of RgbImage Class
*      \date Documented on 10/08/2010
*
*/



/*! RgbImage constructor
	*	\param numRows number of rows
	*   \param numCols number of cols
	*/
inline RgbImage::RgbImage( int numRows, int numCols )
{
	NumRows = numRows;
	NumCols = numCols;
	ImagePtr = new unsigned char[NumRows*GetNumBytesPerRow()];
	if ( !ImagePtr ) {
		fprintf(stderr, "Unable to allocate memory for %ld x %ld bitmap.\n", 
				NumRows, NumCols);
		Reset();
		ErrorCode = MemoryError;
	}
	// Zero out the image
	unsigned char* c = ImagePtr;
	int rowLen = GetNumBytesPerRow();
	for ( int i=0; i<NumRows; i++ ) 
	{
		for ( int j=0; j<rowLen; j++ ) 
		{
			*(c++) = 0;
		}
	}
}

/* ********************************************************************
 *  LoadBmpFile
 *  Read into memory an RGB image from an uncompressed BMP file.
 *  Return true for success, false for failure.  Error code is available
 *     with a separate call.
 *  Author: Sam Buss December 2001.
 **********************************************************************/
/*!
	 *	LoadBmpFile 
	 * \param filename input file name
	 * \return whether the loading is successful
	 */
inline bool RgbImage::LoadBmpFile( const char* filename ) 
{  
	Reset();
	FILE* infile = fopen( filename, "rb" );		// Open for reading binary data
	if ( !infile ) {
		fprintf(stderr, "Unable to open file: %s\n", filename);
		ErrorCode = OpenError;
		return false;
	}

	bool fileFormatOK = false;
	int bChar = fgetc( infile );
	int mChar = fgetc( infile );
	if ( bChar=='B' && mChar=='M' ) {			// If starts with "BM" for "BitMap"
		skipChars( infile, 4+2+2+4+4 );			// Skip 4 fields we don't care about
		NumCols = readLong( infile );
		NumRows = readLong( infile );
		skipChars( infile, 2 );					// Skip one field
		int bitsPerPixel = readShort( infile );
		skipChars( infile, 4+4+4+4+4+4 );		// Skip 6 more fields

		if ( NumCols>0 && NumCols<=100000 && NumRows>0 && NumRows<=100000  
			&& bitsPerPixel==24 && !feof(infile) ) {
			fileFormatOK = true;
		}
	}
	if ( !fileFormatOK ) {
		Reset();
		ErrorCode = FileFormatError;
		fprintf(stderr, "Not a valid 24-bit bitmap file: %s.\n", filename);
		fclose ( infile );
		return false;
	}

	// Allocate memory
	ImagePtr = new unsigned char[NumRows*GetNumBytesPerRow()];
	if ( !ImagePtr ) {
		fprintf(stderr, "Unable to allocate memory for %ld x %ld bitmap: %s.\n", 
				NumRows, NumCols, filename);
		Reset();
		ErrorCode = MemoryError;
		fclose ( infile );
		return false;
	}

	unsigned char* cPtr = ImagePtr;
	for ( int i=0; i<NumRows; i++ ) {
		int j;
		for ( j=0; j<NumCols; j++ ) {
			*(cPtr+2) = fgetc( infile );	// Blue color value
			*(cPtr+1) = fgetc( infile );	// Green color value
			*cPtr = fgetc( infile );		// Red color value
			cPtr += 3;
		}
		int k=3*NumCols;					// Num bytes already read
		for ( ; k<GetNumBytesPerRow(); k++ ) {
			fgetc( infile );				// Read and ignore padding;
			*(cPtr++) = 0;
		}
	}
	if ( feof( infile ) ) {
		fprintf( stderr, "Premature end of file: %s.\n", filename );
		Reset();
		ErrorCode = ReadError;
		fclose ( infile );
		return false;
	}
	fclose( infile );	// Close the file
	return true;
}
/*!	read a short integer form a file
	\param infile input file
	\return short integer just read
	*/
inline short RgbImage::readShort( FILE* infile )
{
	// read a 16 bit integer
	unsigned char lowByte, hiByte;
	lowByte = fgetc(infile);			// Read the low order byte (little endian form)
	hiByte = fgetc(infile);			// Read the high order byte

	// Pack together
	short ret = hiByte;
	ret <<= 8;
	ret |= lowByte;
	return ret;
}
/*!	read a long integer form a file
	\param infile input file
	\return long integer just read
	*/
inline long RgbImage::readLong( FILE* infile )
{  
	// Read in 32 bit integer
	unsigned char byte0, byte1, byte2, byte3;
	byte0 = fgetc(infile);			// Read bytes, low order to high order
	byte1 = fgetc(infile);
	byte2 = fgetc(infile);
	byte3 = fgetc(infile);

	// Pack together
	long ret = byte3;
	ret <<= 8;
	ret |= byte2;
	ret <<= 8;
	ret |= byte1;
	ret <<= 8;
	ret |= byte0;
	return ret;
}
/*!	skip several chars in the file
	\param numChars the number of chars to be skipped
	*/
	
inline void RgbImage::skipChars( FILE* infile, int numChars )
{
	for ( int i=0; i<numChars; i++ ) {
		fgetc( infile );
	}
}

/* ********************************************************************
 *  WriteBmpFile
 *  Write an RGB image to an uncompressed BMP file.
 *  Return true for success, false for failure.  Error code is available
 *     with a separate call.
 *  Author: Sam Buss, January 2003.
 **********************************************************************/
/*!
	 *	WriteBmpFile 
	 * \param filename output file name
	 * \return whether the writing is successful
	 */
inline bool RgbImage::WriteBmpFile( const char* filename )
{
	FILE* outfile = fopen( filename, "wb" );		// Open for reading binary data
	if ( !outfile ) {
		fprintf(stderr, "Unable to open file: %s\n", filename);
		ErrorCode = OpenError;
		return false;
	}

	fputc('B',outfile);
	fputc('M',outfile);
	int rowLen = GetNumBytesPerRow();
	writeLong( 40+14+NumRows*rowLen, outfile );	// Length of file
	writeShort( 0, outfile );					// Reserved for future use
	writeShort( 0, outfile );
	writeLong( 40+14, outfile );				// Offset to pixel data
	writeLong( 40, outfile );					// header length
	writeLong( NumCols, outfile );				// width in pixels
	writeLong( NumRows, outfile );				// height in pixels (pos for bottom up)
	writeShort( 1, outfile );		// number of planes
	writeShort( 24, outfile );		// bits per pixel
	writeLong( 0, outfile );		// no compression
	writeLong( 0, outfile );		// not used if no compression
	writeLong( 0, outfile );		// Pixels per meter
	writeLong( 0, outfile );		// Pixels per meter
	writeLong( 0, outfile );		// unused for 24 bits/pixel
	writeLong( 0, outfile );		// unused for 24 bits/pixel

	// Now write out the pixel data:
	unsigned char* cPtr = ImagePtr;
	for ( int i=0; i<NumRows; i++ ) {
		// Write out i-th row's data
		int j;
		for ( j=0; j<NumCols; j++ ) {
			fputc( *(cPtr+2), outfile);		// Blue color value
			fputc( *(cPtr+1), outfile);		// Blue color value
			fputc( *(cPtr+0), outfile);		// Blue color value
			cPtr+=3;
		}
		// Pad row to word boundary
		int k=3*NumCols;					// Num bytes already read
		for ( ; k<GetNumBytesPerRow(); k++ ) {
			fputc( 0, outfile );				// Read and ignore padding;
			cPtr++;
		}
	}

	fclose( outfile );	// Close the file
	return true;
}
/*!	write a long integer to the file
	\param data the long integter to write
	\param outfile output file
	*/
inline void RgbImage::writeLong( long data, FILE* outfile )
{  
	// Read in 32 bit integer
	unsigned char byte0, byte1, byte2, byte3;
	byte0 = (unsigned char)(data&0x000000ff);		// Write bytes, low order to high order
	byte1 = (unsigned char)((data>>8)&0x000000ff);
	byte2 = (unsigned char)((data>>16)&0x000000ff);
	byte3 = (unsigned char)((data>>24)&0x000000ff);

	fputc( byte0, outfile );
	fputc( byte1, outfile );
	fputc( byte2, outfile );
	fputc( byte3, outfile );
}
/*!	write a short integer to the file
	\param data the short integter to write
	\param outfile output file
	*/
inline void RgbImage::writeShort( short data, FILE* outfile )
{  
	// Read in 32 bit integer
	unsigned char byte0, byte1;
	byte0 = data&0x000000ff;		// Write bytes, low order to high order
	byte1 = (data>>8)&0x000000ff;

	fputc( byte0, outfile );
	fputc( byte1, outfile );
}


/*********************************************************************
 * SetRgbPixel routines allow changing the contents of the RgbImage. *
 *********************************************************************/
/*! set the pixel at (row, col)
	 * \param row row position
	 * \param col column position
	 * \param red   the red color double value
	 * \param green the green color double value
	 * \param blue  the blue  color double value
	 */
inline void RgbImage::SetRgbPixelf( long row, long col, double red, double green, double blue )
{
	SetRgbPixelc( row, col, doubleToUnsignedChar(red), 
							doubleToUnsignedChar(green),
							doubleToUnsignedChar(blue) );
}
/*! set the pixel at (row, col)
	 * \param row row position
	 * \param col column position
	 * \param red   the red color unsigned char
	 * \param green the green color unsigned char
	 * \param blue  the blue  color unsigned char
	 */
inline void RgbImage::SetRgbPixelc( long row, long col,
				   unsigned char red, unsigned char green, unsigned char blue )
{
	assert ( row<NumRows && col<NumCols );
	unsigned char* thePixel = GetRgbPixel( row, col );
	*(thePixel++) = red;
	*(thePixel++) = green;
	*(thePixel) = blue;
}

/*!	convert a double to unsigned char
	* \param x input double number
	* \return output unsigned char
	*/

inline unsigned char RgbImage::doubleToUnsignedChar( double x )
{
	if ( x>=1.0 ) {
		return (unsigned char)255;
	}
	else if ( x<=0.0 ) {
		return (unsigned char)0;
	}
	else {
		return (unsigned char)(x*255.0);		// Rounds down
	}
}
// Bitmap file format  (24 bit/pixel form)		BITMAPFILEHEADER
// Header (14 bytes)
//	 2 bytes: "BM"
//   4 bytes: long int, file size
//   4 bytes: reserved (actually 2 bytes twice)
//   4 bytes: long int, offset to raster data
// Info header (40 bytes)						BITMAPINFOHEADER
//   4 bytes: long int, size of info header (=40)
//	 4 bytes: long int, bitmap width in pixels
//   4 bytes: long int, bitmap height in pixels
//   2 bytes: short int, number of planes (=1)
//   2 bytes: short int, bits per pixel
//   4 bytes: long int, type of compression (not applicable to 24 bits/pixel)
//   4 bytes: long int, image size (not used unless compression is used)
//   4 bytes: long int, x pixels per meter
//   4 bytes: long int, y pixels per meter
//   4 bytes: colors used (not applicable to 24 bit color)
//   4 bytes: colors important (not applicable to 24 bit color)
// "long int" really means "unsigned long int"
// Pixel data: 3 bytes per pixel: RGB values (in reverse order).
//	Rows padded to multiples of four.


#ifndef RGBIMAGE_DONT_USE_OPENGL
/*! Load the Bmp Image from OpenGL buffer
*/
inline bool RgbImage::LoadFromOpenglBuffer()					// Load the bitmap from the current OpenGL buffer
{
	int viewportData[4];
	glGetIntegerv( GL_VIEWPORT, viewportData );
	int& vWidth = viewportData[2];
	int& vHeight = viewportData[3];
	
	if ( ImagePtr==0 ) { // If no memory allocated
		NumRows = vHeight;
		NumCols = vWidth;
		ImagePtr = new unsigned char[NumRows*GetNumBytesPerRow()];
		if ( !ImagePtr ) {
			fprintf(stderr, "Unable to allocate memory for %ld x %ld buffer.\n", 
					NumRows, NumCols);
			Reset();
			ErrorCode = MemoryError;
			return false;
		}
	}
	assert ( vWidth>=NumCols && vHeight>=NumRows );
	int oldGlRowLen;
	if ( vWidth>=NumCols ) {
		glGetIntegerv( GL_UNPACK_ROW_LENGTH, &oldGlRowLen );
		glPixelStorei( GL_UNPACK_ROW_LENGTH, NumCols );
	}
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	// Get the frame buffer data.
	glReadPixels( 0, 0, NumCols, NumRows, GL_RGB, GL_UNSIGNED_BYTE, ImagePtr);

	// Restore the row length in glPixelStorei  (really ought to restore alignment too).
	if ( vWidth>=NumCols ) {
		glPixelStorei( GL_UNPACK_ROW_LENGTH, oldGlRowLen );
	}	
	return true;
}

#endif   // RGBIMAGE_DONT_USE_OPENGL

#endif // RGBIMAGE_H
