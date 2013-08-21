#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LookUpTable/LookUpTable.h"
#include "Frame/Frame.h"
#include "bmpfile/bmpfile.h"
#include "FFT/fft.h"

#ifdef WIN32 // Mkdir
#include <direct.h>
#else
#include <sys/stat.h>
#define _mkdir(x) mkdir(x, 01777)
#endif


#define CURVE_WIDTH 1024
#define CURVE_HEIGHT 512

static LookUpTable * sqrt_table = NULL;

typedef struct
{
	char *name;
	char *path;
	char *output_folder;
	long int size;
	
	bmpfile_t *bmp;
	bmpfile_t *distance_curve;
	bmpfile_t *intensity;
	FILE *handler;
	Frame *frame;

} BinaryData;

long int get_filesize (FILE *f)
{
	fseek (f, 0, SEEK_END);
	long int size = ftell (f);
	rewind (f);
	return size;
}

COMPLEX ** frame_to_complex (Frame *frame)
{
	COMPLEX **c = malloc (sizeof (COMPLEX *) * frame->size);

	for (int i = 0; i < frame->size; i++)
	{
		c[i] = malloc (sizeof (COMPLEX) * frame->size);
	}

	for (int y = 0; y < frame->size; y++)
	{
		for (int x = 0; x < frame->size; x++)
		{
			c[y][x].imag = 0;
			c[y][x].real = frame->data[y][x];
		}
	}

	return c;
}

bmpfile_t * complex_to_bmp (COMPLEX **c, int size)
{
	bmpfile_t *bmp = bmp_create (size, size, 8);

	for (int y = 0; y < size; y++)
	{
		for (int x = 0; x < size; x++)
		{
			unsigned char r = (c[y][x].real * 256.0 > 255.0) ? 
				255 : (unsigned char)(c[y][x].real * 256.0);
			unsigned char b = (c[y][x].imag * 256.0 > 255.0) ? 
				255 : (unsigned char)(c[y][x].imag * 256.0);
			unsigned char g = ((c[y][x].real * 128.0) + (c[y][x].imag * 128.0) > 255) ?
				255 : (unsigned char)(c[y][x].real * 128.0) + (c[y][x].imag * 128.0);

			rgb_pixel_t pixel = {r, g, b, 0};
			bmp_set_pixel (bmp, x, y, pixel);
		}
	}

	return bmp;
}

char *get_filename (char *pathname)
{
	char *res = strrchr (pathname, '/');

	if (res != NULL)
		return res + 1;
	else
		return pathname;
}


char *remove_ext (char* mystr)
{
	char dot = '.';
	char sep = '/';
	// (paxdiablo) http://stackoverflow.com/questions/2736753/how-to-remove-extension-from-file-name

    char *retstr, *lastdot, *lastsep;

    // Error checks and allocate string.
    if (mystr == NULL)
        return NULL;
    if ((retstr = malloc (strlen (mystr) + 1)) == NULL)
        return NULL;

    // Make a copy and find the relevant characters.

    strcpy (retstr, mystr);
    lastdot = strrchr (retstr, dot);
    lastsep = (sep == 0) ? NULL : strrchr (retstr, sep);

    // If it has an extension separator.

    if (lastdot != NULL) {
        // and it's before the extenstion separator.

        if (lastsep != NULL) {
            if (lastsep < lastdot) {
                // then remove it.

                *lastdot = '\0';
            }
        } else {
            // Has extension separator with no path separator.

            *lastdot = '\0';
        }
    }

    // Return the modified string.

    return retstr;
}

double get_euclidian_distance (unsigned char src_x, unsigned char src_y, unsigned char dest_x, unsigned char dest_y)
{
	return lookup_table_get_safe (sqrt_table,
	    pow (dest_x - src_x, 2)
	+   pow (dest_y - src_y, 2)
	);
}

void binary_save_bmp (BinaryData *binary, unsigned long int start_offset, unsigned long int cur_offset)
{
	char bmp_filename[1024];
	static unsigned int frame_id = 0;

	sprintf (bmp_filename, "%s/%s/%d_%s_0x%lx-0x%lx.bmp", 
		binary->output_folder, binary->path, frame_id++, binary->name, start_offset, cur_offset);

	// Output to bmp
	bmp_save (binary->bmp, bmp_filename);
	
	// Cleaning
	bmp_destroy (binary->bmp);
	binary->bmp = bmp_create (256, 256, 8);
}

void binary_save_frame (BinaryData *binary, unsigned long int start_offset, unsigned long int cur_offset)
{
	char bmp_filename[1024];
	static unsigned int frame_id = 0;

	sprintf (bmp_filename, "%s/%s/%d_%s_FFT_0x%lx-0x%lx.bmp", 
		binary->output_folder, binary->path, frame_id++, binary->name, start_offset, cur_offset);

	// Convert frame to complex matrix and compute FFT
	COMPLEX **c = frame_to_complex (binary->frame);
	FFT2D (c, binary->frame->size, binary->frame->size, 1);
	
	// Output to bmp
	bmpfile_t *bmp = complex_to_bmp (c, binary->frame->size);
	bmp_save (bmp, bmp_filename);
	
	// Cleaning
	bmp_destroy (bmp);
	frame_reset (binary->frame);
}

void binary_save (BinaryData *binary, unsigned long int start_offset, unsigned long int cur_offset)
{
	binary_save_bmp (binary, start_offset, cur_offset);
	binary_save_frame (binary, start_offset, cur_offset);
}


void bmp_save_distance_curve (BinaryData *binary)
{
	char bmp_filename[1024];

	sprintf (bmp_filename, "%s/curves/%s_distance_curve.bmp", binary->output_folder, binary->name);
	bmp_save (binary->distance_curve, bmp_filename);
	bmp_destroy (binary->distance_curve);
}

double compute_moyenne (double sum, unsigned int count)
{
	return sum / count;
}

double compute_variance (double sq_sum, double moyenne, int n)
{
	return (sq_sum / n) - pow (moyenne, 2);
}

double compute_ecart_type (double moyenne, double sq_sum, int n)
{
    double variance = compute_variance (sq_sum, moyenne, n);
    return lookup_table_get_safe (sqrt_table, variance);
}

double compute_rsd (double ecart_type, double moyenne)
{
	if (moyenne == 0.0)
		return 0.0;

	return ecart_type / moyenne;
}

void bresenham (BinaryData *binary, int x1, int y1, int x2, int y2, rgb_pixel_t pixel)
{
	int dx = x2 - x1, dy = y2 - y1, i, e;
	int incx = 1, incy = 1, inc1, inc2;
	int x, y;

	if (dx < 0)  dx = -dx;
	if (dy < 0)  dy = -dy;
	if (x2 < x1) incx = -1;
	if (y2 < y1) incy = -1;

	x = x1;
	y = y1;

	if (dx > dy)
	{
		bmp_set_pixel (binary->distance_curve, x, y, pixel);
		e = 2 * dy - dx;
		inc1 = 2 * (dy -dx);
		inc2 = 2 * dy;

		for (i = 0; i < dx; i++)
		{
			if (e >= 0)
			{
				y += incy;
				e += inc1;
			}

			else
				e += inc2;

			x += incx;
			bmp_set_pixel (binary->distance_curve, x, y, pixel);
		}
	}

	else
	{
		bmp_set_pixel (binary->distance_curve, x, y, pixel);
		e = 2 * dx - dy;
		inc1 = 2 * (dx - dy);
		inc2 = 2 * dx;

		for (i = 0; i < dy; i++)
		{
			if (e >= 0)
			{
				x += incx;
				e += inc1;
			}

			else
				e += inc2;

			y += incy;
			bmp_set_pixel (binary->distance_curve, x, y, pixel);
		}
	}
}

void binary_draw_pixel (BinaryData *binary, unsigned char pixel_pos[2], rgb_pixel_t pixel)
{
	frame_inc (binary->frame, pixel_pos[0], pixel_pos[1]);
	bmp_set_pixel (binary->bmp, pixel_pos[0], pixel_pos[1], pixel);
}

void analyze_2d (BinaryData *binary)
{
	printf("Analyze 2D started (%s : %ld bytes)\n", binary->name, binary->size);

	// Pixel
	unsigned char pixel_pos[2];
	unsigned char pixel_last_pos[2];
	unsigned char pixel_color;

	// File Offsets
	long int cur_offset   = 0;
	long int start_offset = 0;
	int first_iteration   = 1;
	int first_curve_point = 1;
	unsigned int c;

	// Frame related
	unsigned int points_in_frame = 0;
	unsigned int points_far = 0;
	double frame_distance = 0.0;
	double frame_distance_square = 0.0;
	double last_distance = -1.0;
	// Distance curve
	double last_curve_progress   = 0.0;
	double last_curve_ecart_type = -1.0;
	double last_curve_moyenne    = -1.0;
	double last_curve_rsd        = -1.0;
	unsigned int points_in_curve = 0;
	double curve_distance_square = 0.0;
	double curve_distance        = 0.0;

	// Colors
	rgb_pixel_t red   = {0,   0,   255, 0};
	rgb_pixel_t blue  = {255, 128, 128, 0};
	rgb_pixel_t green = {128, 255, 128, 0};
	rgb_pixel_t black = {0,   0,   0,   0};

	printf ("     Offset       |    Distance\n");
	pixel_pos[1] = fgetc (binary->handler);
	while ((c = fgetc (binary->handler)) != EOF)
	{
		// Cortesi algorithm
		pixel_pos[0] = pixel_pos[1];
		pixel_pos[1] = c;
		points_in_frame++;
		points_in_curve++;

		// Color computation
		double progress = ((double) cur_offset / binary->size);
		pixel_color = progress * 256;

		// Draw in bitmap
		rgb_pixel_t pixel = {255 - pixel_color, 128 - (pixel_color / 2.0), pixel_color, 0};
		binary_draw_pixel (binary, pixel_pos, pixel);

		// Get the distance from the last point
		if (first_iteration)
		{
			pixel_last_pos[0] = pixel_pos[0];
			pixel_last_pos[1] = pixel_pos[1];
		}

		double distance = get_euclidian_distance (pixel_pos[0], pixel_pos[1], pixel_last_pos[0], pixel_last_pos[1]);
		frame_distance += distance;
		frame_distance_square += pow (distance, 2);
		curve_distance += distance;
		curve_distance_square += pow (distance, 2);

		if (first_iteration)
			last_distance = distance;

		double moyenne        = compute_moyenne (frame_distance, points_in_frame);
		double ecart_type     = compute_ecart_type (moyenne, frame_distance_square, points_in_frame);
		double last_moyenne   = compute_moyenne (distance + last_distance, 2);

		// Draw frames
		if (last_moyenne > ecart_type * (moyenne / last_moyenne))
		{
			long int maxsize = binary->size / 10;
			if (maxsize > 256 * 256)
				maxsize = 256 * 256;

			if (points_far++ > maxsize && frame_distance > 800000.0)
			{
				printf ("%.8lx-%.8lx | %f\n", start_offset, cur_offset, frame_distance);
				binary_save (binary, start_offset, cur_offset);

				bresenham (binary,
					 (int)(progress * CURVE_WIDTH), CURVE_HEIGHT,
					 (int)(progress * CURVE_WIDTH), 0,
					black
				);

				points_in_frame = 0;
				frame_distance  = 0.0;
				frame_distance_square = 0.0;
				points_far = 0;
				start_offset = cur_offset;
			}
		}

		// Draw distance curve
		if ((int)(CURVE_WIDTH * progress) != (int)(CURVE_WIDTH * last_curve_progress))
		{
			double curve_moyenne    = compute_moyenne (curve_distance, points_in_curve);
			double curve_ecart_type = compute_ecart_type (curve_moyenne, curve_distance_square, points_in_curve);
			double curve_rsd        = compute_rsd (curve_ecart_type, curve_moyenne) * (CURVE_HEIGHT / 2.0);

			if (first_curve_point)
			{
				last_curve_ecart_type = curve_ecart_type;
				last_curve_moyenne    = curve_moyenne;
				last_curve_rsd        = curve_rsd;
			}

			bresenham (binary,
				 (int)(last_curve_progress * CURVE_WIDTH), (int)(CURVE_HEIGHT - last_curve_ecart_type - 1),
				 (int)(progress            * CURVE_WIDTH), (int)(CURVE_HEIGHT - curve_ecart_type - 1),
				red
			);

			bresenham (binary,
				 (int)(last_curve_progress * CURVE_WIDTH), (int)(CURVE_HEIGHT - last_curve_moyenne - 1),
				 (int)(progress            * CURVE_WIDTH), (int)(CURVE_HEIGHT - curve_moyenne - 1),
				blue
			);

			bresenham (binary,
				 (int)(last_curve_progress * CURVE_WIDTH), (int)(CURVE_HEIGHT - last_curve_rsd - 1),
				 (int)(progress            * CURVE_WIDTH), (int)(CURVE_HEIGHT - curve_rsd - 1),
				green
			);

			last_curve_ecart_type = curve_ecart_type;
			last_curve_moyenne    = curve_moyenne;
			last_curve_rsd        = curve_rsd;
			last_curve_progress   = progress;

			points_in_curve = 0;
			curve_distance  = 0.0;
			curve_distance_square = 0.0;
			first_curve_point = 0;
		}

		// Final computation before next iteration
		pixel_last_pos[0] = pixel_pos[0];
		pixel_last_pos[1] = pixel_pos[1];
		cur_offset++;
		last_distance = distance;
		first_iteration = 0;
	}

	if (points_in_frame > 0)
		binary_save (binary, start_offset, cur_offset);

	bmp_save_distance_curve (binary);
}

int main (int argc, char **argv)
{
	if (argc < 2)
	{
		printf ("Usage : %s <binary> [output folder]\n", argv[0]);
		return 1;
	}

	FILE *bin = fopen (argv[1], "rb");
	if (!bin)
	{
		printf ("File %s cannot be opened.\n", argv[1]);
		return 2;
	}
	
	char *output_folder = (argc >= 3) ? argv[2] : "./";

	// LookUp table initialization
	sqrt_table = lookup_table_new (-0.1, (256 * 256) * 2, 0.1, sqrt);

	// Binary file related
	BinaryData binary = {
		.name = get_filename (argv[1]),
		.path = remove_ext (get_filename (argv[1])),
		.size = get_filesize (bin),
		.bmp  = bmp_create (256, 256, 8),
		.distance_curve = bmp_create (CURVE_WIDTH, CURVE_HEIGHT, 8),
		.intensity = bmp_create (256, 256, 8),
		.handler = bin,
		.frame = frame_new (256),
		.output_folder = output_folder
	};

	// Prepare folders
	char tmp[1024];
	_mkdir (binary.output_folder);
	sprintf(tmp, "%s/%s", binary.output_folder, binary.path);
	_mkdir (tmp);
	sprintf(tmp, "%s/curves", binary.output_folder);
	_mkdir (tmp);


	// Analyze
	analyze_2d (&binary);

	return 0;
}
