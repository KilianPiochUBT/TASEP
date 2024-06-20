#include "graphlib.h"

const pixel_t white = { 255, 255, 255 };
const pixel_t black = { 0, 0, 0 };
const pixel_t red = { 255, 0, 0 };
const pixel_t green = { 0, 255, 0 };
const pixel_t blue = { 0, 0, 255 };

screen_t *create_screen(int rx, int ry)
{
  assert(rx >= 0);
  assert(ry >= 0);
  
  screen_t *screen = (screen_t *) malloc(sizeof(screen_t));
  assert(screen != NULL);

  screen->data = (pixel_t *) malloc(rx * ry * sizeof(pixel_t));
  assert(screen->data != NULL);

  screen->rx =rx;
  screen->ry =ry;

  return screen;
} 
  
void destroy_screen(screen_t *screen)
{
  free(screen->data);
  free(screen);
}

/* Funktion zum Speichern des Inhalts des Bildspeichers als
   PPM-Datei */

void write_ppm(screen_t *screen, const char *filename)
{
  FILE *fp;
  int z;
  pixel_t *p;

  assert(screen != NULL);
  assert(screen->data != NULL);

  z = screen->rx * screen->ry;
  p = screen->data;

  if (filename == NULL)
    fp = stdout;
  else
  {
    fp = fopen(filename, "w");
    assert(fp != NULL);
  }
      
  fprintf(fp, "P6\n%d %d\n255\n", screen->rx, screen->ry);

  while (z--)
    fwrite(p++, 1, 3, fp);

  if (fp != stdout)
    fclose(fp);
}

/* Zeichenfunktionen */

/* F�lle Bildinhalt mit Farbe col */

void clear_screen(screen_t *screen)
{
  fill_screen(screen, white);
}

void fill_screen(screen_t *screen, const pixel_t col)
{
  assert(screen != NULL);
  assert(screen->data != NULL);

  pixel_t *p = screen->data + screen->rx * screen->ry;

  while (screen->data != p--)
    *p = col;
}

/* Setze Farbe des Pixels an Position (x, y) auf col */

void set_pixel(screen_t *screen, int x, int y, pixel_t col)
{
  assert(screen != NULL);
  if (x < 0 || x >= screen->rx || y < 0 || y >= screen->ry)
  {
    //printf("Access to screen[%d][%d]\n", x, y);
  }
  else
  {
    assert(screen->data != NULL);
    screen->data[y * screen->rx + x] = col;
  }
}
