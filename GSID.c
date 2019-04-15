#include "GSID.h"
#include "time.h"

double E = 2.718281828459;
double SIGMA;
int MAX_PATH_LENGTH = 200;
int SAMPLES = 20;
int regionSize = 9;  // MUST BE ODD
Image *img;
double *rgbIm;
Image *normalMap;
double *nmap;
Image *matMap;
double *mmap;
Image *shadowMap;
double *smap;
Image *depthMap;
double *dmap;
Image *colourMap;
double *cmap;
int DEBUGLOOP = 0;
double THRESHOLD;
double *allProbs;
Image *pOut;
double *pathOut;
int mid;

double *hits, *total;

double rgbWeight = 0;
double normalWeight = 0;
double diffWeight = 0;
double reflWeight = 0;
double refrWeight = 0;
double shadowWeight = 0;
double depthWeight = 0;
double colWeight = 0;

// Gets one dimensional index from x and y
// inline int ((int x) + ( int y) { return (x + (y * img->sx)) * img->sx); }

double probWeighted(int idx, int xj) {
  // image colour comparison
  double xi_to_xj = 0;
  double dr = rgbIm[idx + 0] - rgbIm[xj + 0];
  double dg = rgbIm[idx + 1] - rgbIm[xj + 1];
  double db = rgbIm[idx + 2] - rgbIm[xj + 2];
  xi_to_xj += rgbWeight * (dr * dr + dg * dg + db * db);

  // dot product comparison
  double n1x = 2 * nmap[idx] - 1;
  double n1y = 2 * nmap[idx + 1] - 1;
  double n1z = 2 * nmap[idx + 2] - 1;
  double n2x = 2 * nmap[xj] - 1;
  double n2y = 2 * nmap[xj + 1] - 1;
  double n2z = 2 * nmap[xj + 2] - 1;
  double dotProd = n1x * n2x + n1y * n2y + n1z * n2z;
  dotProd = (dotProd + 1) / 2;
  dotProd = 1 - dotProd;
  xi_to_xj += normalWeight * dotProd;
  // added plus TOL to the differences because if any of those features are
  // running solo a lot of the time all 8 directions have prob 0

  // material comparison
  double diffComp = fabs(mmap[idx] - mmap[xj]);
  xi_to_xj += diffWeight * diffComp;
  double reflComp = fabs(mmap[idx + 1] - mmap[xj + 1]);
  xi_to_xj += reflWeight * reflComp;
  double refrComp = fabs(mmap[idx + 2] - mmap[xj + 2]);
  xi_to_xj += refrWeight * refrComp;

  // depth comparison
  double depthComp = fabs(dmap[idx] - dmap[xj]);
  xi_to_xj += depthWeight * depthComp;

  // shadow comparison
  // TODO maybe check surrounding pixels and take an average since shadows and
  // caustics are noisy in path tracing
  double shadowComp = fabs(smap[idx] - smap[xj]);
  xi_to_xj += shadowWeight * shadowComp;

  dr = fabs(cmap[idx + 0] - cmap[xj + 0]);
  dg = fabs(cmap[idx + 1] - cmap[xj + 1]);
  db = fabs(cmap[idx + 2] - cmap[xj + 2]);
  xi_to_xj += colWeight * (dr * dr + dg * dg + db * db);

  return -xi_to_xj * SIGMA;
}

// Note that the probability return here is
// NOT normalized. Normalize all the probs.
// for the neighbours of a pixel together.
double probability(int x0, int xj, int xj1, int dir, double *rootProbs, int x,
                   int y, int i, int j) {
  int c0 = 3 * x0;
  int cj = 3 * xj;
  int cj1 = 3 * xj1;
  double term_1_exponent;
  // if xj1 within n sized grid of x0 get stored value other wise compute it
  x = x - i + mid;
  y = y - j + mid;
  if (x < 0 || y < 0 || x >= regionSize || y >= regionSize) {
    term_1_exponent = probWeighted(c0, cj1);
  } else {
    hits[x0] = hits[x0] + 1;
    term_1_exponent = rootProbs[x + y * regionSize];
  }
  total[x0] = total[x0] + 1;

  // double term_1_exponent = probWeighted(c0, cj1);
  double term_2_exponent = allProbs[xj + dir];  //-xj_to_xj1 * SIGMA;

  // Working with logs... finally.
  return term_1_exponent + term_2_exponent;
}
// uses the above probability function to find the probability of randomly
// walking in potentially all 4 directions, then normalizes them
void generateProbs(double probs[8], int x, int y, int root, double *rootProbs,
                   int i, int j) {
  // Initialize them
  probs[0] = probs[1] = probs[2] = probs[3] = probs[4] = probs[5] = probs[6] =
      probs[7] = 100;

  int CUR = ((x) + (y)*img->sx);

  if (x + 1 < img->sx) {
    int RIGHT = ((x + 1) + (y)*img->sx);
    probs[0] = probability(root, CUR, RIGHT, 0, rootProbs, x + 1, y, i, j);
  }
  if (x - 1 >= 0) {
    int LEFT = ((x - 1) + (y)*img->sx);
    probs[1] = probability(root, CUR, LEFT, 1, rootProbs, x - 1, y, i, j);
  }
  if (y + 1 < img->sy) {
    int BOTTOM = ((x) + (y + 1) * img->sx);
    probs[2] = probability(root, CUR, BOTTOM, 2, rootProbs, x, y + 1, i, j);
  }
  if (y - 1 >= 0) {
    int TOP = ((x) + (y - 1) * img->sx);
    probs[3] = probability(root, CUR, TOP, 3, rootProbs, x, y - 1, i, j);
  }
  if (x + 1 < img->sx && y + 1 < img->sy) {
    int BOTTOM_RIGHT = ((x + 1) + (y + 1) * img->sx);
    probs[4] =
        probability(root, CUR, BOTTOM_RIGHT, 4, rootProbs, x + 1, y + 1, i, j);
  }
  if (x - 1 >= 0 && y + 1 < img->sy) {
    int BOTTOM_LEFT = ((x - 1) + (y + 1) * img->sx);
    probs[5] =
        probability(root, CUR, BOTTOM_LEFT, 5, rootProbs, x - 1, y + 1, i, j);
  }
  if (x + 1 < img->sx && y - 1 >= 0) {
    int TOP_RIGHT = ((x + 1) + (y - 1) * img->sx);
    probs[6] =
        probability(root, CUR, TOP_RIGHT, 6, rootProbs, x + 1, y - 1, i, j);
  }
  if (x - 1 >= 0 && y - 1 >= 0) {
    int TOP_LEFT = ((x - 1) + (y - 1) * img->sx);
    probs[7] =
        probability(root, CUR, TOP_LEFT, 7, rootProbs, x - 1, y - 1, i, j);
  }

  // generate normalization constant and normalize the probabilities
  double sum = 0.0;
  if (probs[0] != 100) sum += exp(probs[0]);
  if (probs[1] != 100) sum += exp(probs[1]);
  if (probs[2] != 100) sum += exp(probs[2]);
  if (probs[3] != 100) sum += exp(probs[3]);
  if (probs[4] != 100) sum += exp(probs[4]);
  if (probs[5] != 100) sum += exp(probs[5]);
  if (probs[6] != 100) sum += exp(probs[6]);
  if (probs[7] != 100) sum += exp(probs[7]);

  // Normalize values... since we're working with logs,
  // we subtract the log of the normalization const.
  double logSum = log(sum);
  if (probs[0] != 100) probs[0] -= logSum;
  if (probs[1] != 100) probs[1] -= logSum;
  if (probs[2] != 100) probs[2] -= logSum;
  if (probs[3] != 100) probs[3] -= logSum;
  if (probs[4] != 100) probs[4] -= logSum;
  if (probs[5] != 100) probs[5] -= logSum;
  if (probs[6] != 100) probs[6] -= logSum;
  if (probs[7] != 100) probs[7] -= logSum;
}

void generateTransitions(double probs[8], int idx) {
  probs[0] = probs[1] = probs[2] = probs[3] = probs[4] = probs[5] = probs[6] =
      probs[7] = 100;
  int x = idx % img->sx;
  int y = idx / img->sx;
  idx = 3 * idx;

  double xi_to_xj;
  if (x + 1 < img->sx) {
    int RIGHT = ((x + 1) + (y)*img->sx) * 3;
    probs[0] = probWeighted(idx, RIGHT);
  }
  if (x - 1 >= 0) {
    int LEFT = ((x - 1) + (y)*img->sx) * 3;
    probs[1] = probWeighted(idx, LEFT);
  }
  if (y + 1 < img->sy) {
    int BOTTOM = ((x) + (y + 1) * img->sx) * 3;
    probs[2] = probWeighted(idx, BOTTOM);
  }
  if (y - 1 >= 0) {
    int TOP = ((x) + (y - 1) * img->sx) * 3;
    probs[3] = probWeighted(idx, TOP);
  }
  if (x + 1 < img->sx && y + 1 < img->sy) {
    int BOTTOM_RIGHT = ((x + 1) + (y + 1) * img->sx) * 3;
    probs[4] = probWeighted(idx, BOTTOM_RIGHT);
  }
  if (x - 1 >= 0 && y + 1 < img->sy) {
    int BOTTOM_LEFT = ((x - 1) + (y + 1) * img->sx) * 3;
    probs[5] = probWeighted(idx, BOTTOM_LEFT);
  }
  if (x + 1 < img->sx && y - 1 >= 0) {
    int TOP_RIGHT = ((x + 1) + (y - 1) * img->sx) * 3;
    probs[6] = probWeighted(idx, TOP_RIGHT);
  }
  if (x - 1 >= 0 && y - 1 >= 0) {
    int TOP_LEFT = ((x - 1) + (y - 1) * img->sx) * 3;
    probs[7] = probWeighted(idx, TOP_LEFT);
  }
}

// Modifies x and y to have the coords of the
// correct neighbour according to 'i'.
//    0 = right
//    1 = left
//    2 = bottom
//    3 = top
void getDir(int *x, int *y, int i) {}

Image *denoise() {
  // Create the denoised image
  Image *new = newImage(img->sx, img->sy);
  double *rgbDenoise = (double *)new->rgbdata;

  // We will store the weights and paths for
  // each of the random walks in these arrays.
  int root;
  int bound = new->sx * new->sy;
  hits = (double *)malloc(sizeof(double) * bound);
  total = (double *)malloc(sizeof(double) * bound);
  for (int i = 0; i < bound; i++) {
    hits[i] = 0;
    total[i] = 0;
  }
  int regionArea = regionSize * regionSize;
  mid = (regionSize - 1) / 2;  /// only for odd regionSize values

  fprintf(stderr, "Denoising row: \n");
#pragma omp parallel for schedule(dynamic, 32) private(root)
  for (root = 0; root < bound; root++) {
    // 1-dimensional index of the current pixel
    int i = root % new->sx;
    int j = root / new->sx;

    // compute region of values for root pixel
    double rootProbs[regionArea];
    for (int b = 0; b < regionArea; b++) {
      int x = b % regionSize;
      int y = b / regionSize;

      x = x + i - mid;
      y = y + j - mid;

      if (x < 0 || y < 0 || x >= img->sx || y >= img->sy) {
        rootProbs[b] = 0;
        continue;
      }

      rootProbs[b] = probWeighted(3 * root, 3 * (x + y * img->sx));
    }

    // SAMPLES amount of random walks per pixel
    double sum_weights = 0;
    double R = 0, G = 0, B = 0;
    for (int m = 0; m < SAMPLES; m++) {
      // Starting pixel coords for the path
      int x = i;
      int y = j;

      double cur_path_prob = 0.0;
      // Max path size is LENGTH
      // prob of Wij for k = 0 should be 1 since it's p(x|x)
      for (int k = 1; k < MAX_PATH_LENGTH; k++) {
        if ((i == 256 && j == 256) || (i == 128 && j == 384) ||
            (i == 1 && j == 1) || (i == 100 && j == 487) ||
            (i == 265 && j == 380)) {
          pathOut[3 * ((x) + (y)*img->sx) + 1] += .1;
          pathOut[3 * ((x) + (y)*img->sx) + 0] =
              pathOut[3 * ((x) + (y)*img->sx) + 1];
          pathOut[3 * ((x) + (y)*img->sx) + 2] = 0;
        }

        // Compute the probabilities for transition
        // to the neighboring pixels and pick one
        double probs[8];
        generateProbs(probs, x, y, root, rootProbs, i, j);
        int dir = randomDir(probs);

        cur_path_prob += probs[dir];

        if (cur_path_prob <= THRESHOLD) {
          break;
        }

        double wgt = exp(cur_path_prob / k);
        sum_weights += wgt;

        // If the probabilities were too small, end path
        if (dir < 0) break;
        // Change current pixel to the neighbour, repeat
        x += (dir == 0 || dir == 4 || dir == 6);
        x -= (dir == 1 || dir == 5 || dir == 7);
        y += (dir == 2 || dir == 4 || dir == 5);
        y -= (dir == 3 || dir == 6 || dir == 7);

        // Store the weight for now, we will accumulate
        // the color after we perform all the walks.
        // if (k == 0) fprintf(stderr, "%f\n", probs[dir]);
        R += rgbIm[3 * ((x) + (y)*img->sx) + 0] * wgt;
        G += rgbIm[3 * ((x) + (y)*img->sx) + 1] * wgt;
        B += rgbIm[3 * ((x) + (y)*img->sx) + 2] * wgt;
      }
    }
    // Set the denoised color estimate for this pixel
    // fprintf(stderr, "Colours: %f, %f, %f\n", R / C, G / C, B / C);
    rgbDenoise[3 * root + 0] = R / sum_weights;
    rgbDenoise[3 * root + 1] = G / sum_weights;
    rgbDenoise[3 * root + 2] = B / sum_weights;

    if (i == 0 && j % 5 == 0) {
      printf("\r%d / %d", j, new->sx);
      fflush(stdout);
    }
  }
  printf("\nDone!\n");

  return new;
}

int main(int argc, char **argv) {
  srand(time(NULL));
  if (argc < 10) {
    fprintf(stderr,
            "Usage: %s input output_name sigma threshold normal_map depth_map "
            "shadow_map material_map colour_map\n",
            argv[0]);
    exit(1);
  }

  SIGMA = atof(argv[3]);
  if (SIGMA < 0 || SIGMA > 1) {
    fprintf(stderr, "Invalid SIGMA. Please enter a value between 0 and 1.\n");
    exit(1);
  }
  // Get threshold and log it so we can do
  // log-based comparisons directly.
  THRESHOLD = atof(argv[4]);
  THRESHOLD = log(THRESHOLD);

  SIGMA = 2 * SIGMA * SIGMA;
  SIGMA = 1 / SIGMA;

  img = readPPMimage(argv[1]);
  normalMap = readPPMimage(argv[5]);
  depthMap = readPPMimage(argv[6]);
  shadowMap = readPPMimage(argv[7]);
  matMap = readPPMimage(argv[8]);
  colourMap = readPPMimage(argv[9]);
  rgbIm = (double *)img->rgbdata;
  nmap = (double *)normalMap->rgbdata;
  dmap = (double *)depthMap->rgbdata;
  cmap = (double *)colourMap->rgbdata;
  mmap = (double *)matMap->rgbdata;
  smap = (double *)shadowMap->rgbdata;
  rgbWeight = 2.0;
  normalWeight = 3.0;
  diffWeight = 1.0;
  reflWeight = 2.0;
  refrWeight = 1.0;
  shadowWeight = 0.5;
  depthWeight = 2.0;
  colWeight = 2.0;

  fprintf(stderr, "Precalculating transitional probabilities...\n");
  int pixels = img->sx * img->sy;
  pOut = newImage(img->sx, img->sy);
  pathOut = (double *)pOut->rgbdata;
  for (int i = 0; i < pixels; i++) {
    pathOut[(3 * i)] = rgbIm[3 * i];
    pathOut[(3 * i) + 1] = rgbIm[(3 * i) + 1];
    pathOut[(3 * i) + 2] = rgbIm[(3 * i) + 2];
  }

  allProbs = (double *)malloc(sizeof(double) * pixels * 8);
  if (!allProbs) exit(0);
  int i;
#pragma omp parallel for schedule(dynamic, 32) private(i)
  for (i = 0; i < pixels; i++) {
    double probs[8];
    int iInd = i * 8;
    generateTransitions(probs, i);
    for (int j = 0; j < 8; j++) {
      *(allProbs + iInd + j) = probs[j];  // xi_to_xj;
    }
  }
  fprintf(stderr, "Done calculating probabilities!\n");

  Image *denoised = denoise();
  double num = 0, den = 0;
  for (int bound = 0; bound < pixels; bound++) {
    num += hits[bound];
    den += total[bound];
  }
  fprintf(stderr, "hit percent: %f", num / den);
  free(hits);
  free(total);
  free(allProbs);

  free(img);
  free(normalMap);
  free(depthMap);
  free(shadowMap);
  free(matMap);
  free(colourMap);

  imageOutput(denoised, argv[2]);
  imageOutput(pOut, "NormalPathGSID");
  free(pOut);
  free(denoised);
  return 0;
}

// Takes in a list of probabilites that
// sum up to 1. Returns an index chosen
// with the probability of that index.
// Returns -1 if probabilites don't add
// up (this happens on early termination)
int randomDir(double *probabilities) {
  double choice = drand48();
  double sum = 0;
  int i;
  for (i = 0; i < 8; i++) {
    if (probabilities[i] != 100) sum += exp(probabilities[i]);
    if (choice < sum) return i;
  }
  return i;
}

// Returns the squared distance between the
// difference of the RGB Values
double distance(Pixel *a, Pixel *b) {
  double dr = (a->r - b->r);
  double dg = (a->g - b->g);
  double db = (a->b - b->b);
  return dr * dr + db * db + dg * dg;
}

// Returns the squared euclidean distance
// between the pixel coords
double distance_euclidean(Pixel *a, Pixel *b) {
  int ax = a->idx % img->sx;
  int ay = a->idx / img->sx;
  int bx = b->idx % img->sx;
  int by = b->idx / img->sx;
  return (ax - bx) * (ax - bx) * 1.0 + (ay - by) * (ay - by) * 1.0;
}
