#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "limace.h"
#include "tai.h"

#define DEBOGAGE

#ifdef DEBOGAGE
#define DEBUG fprintf(stderr,"Fichier %s, ligne %d\n",__FILE__,__LINE__);
#else
#define DEBUG
#endif

#define AFAIRE(ValeurRetour) \
  fprintf(stderr,"--> Fichier %s, ligne %d, corps de la fonction %s à écrire.\n",\
          __FILE__,__LINE__,__func__);\
  return ValeurRetour;  


/*
 * Conversion d'une image couleur en une image de niveaux de gris
 * Entrée : image initiale en couleur
 * Sortie : image de niveaux de gris résultat
 */
Image RGB2Gray(Image Im)
{
  int lignes = ImNbRow(Im);
  int colonnes = ImNbCol(Im);

  Image ImGris = ImAlloc(GrayLevel, lignes, colonnes);

  unsigned char ** I  = ImGetI(ImGris);
  unsigned char ** IR = ImGetR(Im);
  unsigned char ** IG = ImGetG(Im);
  unsigned char ** IB = ImGetB(Im);

  for(int i = 0; i < lignes; i++) {
    for(int j = 0; j < colonnes; j++) {
      I[i][j] = round(0.299 * IR[i][j] + 0.587 * IG[i][j] + 0.114 * IB[i][j]);
    }
  }
  return ImGris;
}


/*
 * Binarisation d'une image de niveaux de gris par seuillage global
 * Entrées : image de niveaux de gris initiale
             seuil (niveau de gris)
 * Sortie : image binaire
 */
Image Binarization(Image Im, unsigned char Threshold)
{
  int lignes = ImNbRow(Im);
  int colonnes = ImNbCol(Im);

  Image ImBin = ImAlloc(BitMap, lignes, colonnes);

  unsigned char ** I = ImGetI(Im);
  unsigned char ** MatBin = ImGetI(ImBin);

  for(int i = 0; i < lignes; i++) {
    for(int j = 0; j < colonnes; j++) {
      MatBin[i][j] = (I[i][j] < Threshold) ? 0 : 1;
    }
  }
  return ImBin;
}


/*
 * Inversion d'une image
 * Entrée : image initiale (binaire, niveaux de gris ou couleur)
 * Sortie : image résultat
 */
Image Inversion(Image Im)
{
  int lignes = ImNbRow(Im);
  int colonnes = ImNbCol(Im);

  ImageType type = ImType(Im);

  Image ImInv = ImAlloc(type, lignes, colonnes);

  if (type != Color) {
    unsigned char ** MatIm = ImGetI(Im);
    unsigned char ** MatInv = ImGetI(ImInv);

    for(int i = 0; i < lignes; i++) {
      for(int j = 0; j < colonnes; j++) {
        if (type == BitMap) {
          MatInv[i][j] = (MatIm[i][j] ? 0 : 1);
        } else {
          MatInv[i][j] = 255 - MatIm[i][j];
        }
      }
    }
  } else {
    unsigned char ** MatRInv = ImGetR(ImInv);
    unsigned char ** MatGInv = ImGetG(ImInv);
    unsigned char ** MatBInv = ImGetB(ImInv);

    unsigned char ** MatRI = ImGetR(Im);
    unsigned char ** MatGI = ImGetG(Im);
    unsigned char ** MatBI = ImGetB(Im);
    for(int i = 0; i < lignes; i++) {
      for(int j = 0; j < colonnes; j++) {
        MatRInv[i][j] = 255-MatRI[i][j];
        MatGInv[i][j] = 255-MatGI[i][j];
        MatBInv[i][j] = 255-MatBI[i][j];
      }
    }
  }
  return ImInv;
}


/*
 * Calcul de l'histogramme d'une image de niveaux de gris
 * Entrée : image initiale (niveaux de gris)
 * Sortie : histogramme (matrice de int 1 x 256)
 */
Matrix Histogram(Image Im)
{
  int lignes = ImNbRow(Im);
  int colonnes = ImNbCol(Im);

  Matrix mat = MatCAlloc(Int,1,256);
  int ** matHist = MatGetInt(mat);

  unsigned char ** matIm = ImGetI(Im);

  for(int i = 0; i < lignes; i++) {
    for(int j = 0; j < colonnes; j++) {
      matHist[0][matIm[i][j]]++;
    }
  }
  return mat;
}


/*
 * Représentation d'un histogramme sous forme d'une image
 * Entrées : histogramme (matrice de int 1 x 256) et nombre de lignes de
 * l'image résultat (une échelle des niveaux de gris de 25 lignes est ajoutée
 * sous l'histogramme)
 * Sortie : image de niveaux de gris résultat
 */
Image Hist2Im(Matrix Hist, int NbLig)
{
	unsigned char **I;
	int *h,i,j,Max=0,NbCol=256,NbLig2=NbLig+25;
	Image Res;

	if (MatType(Hist)!=Int) return NULL;
  NbLig2=NbLig+25;
	Res=ImCAlloc(GrayLevel,NbLig2,NbCol);
  if (Res==NULL) return NULL;
	h=*MatGetInt(Hist);
	for (j=0;j<NbCol;j++)
		if (h[j]>Max) Max=h[j];
	I=ImGetI(Res);
	for (j=0;j<256;j++)
		for (i=NbLig-1;i>=(NbLig-NbLig*h[j]/Max);i--)
		    I[i][j]=255;
  for (j=0;j<256;j++)
    I[NbLig][j]=0;
  for (i=NbLig+1;i<NbLig2;i++)
    for (j=0;j<256;j++)
      I[i][j]=j;
	return Res;
}

/* Fonction Auxiliaire pour le cout d'Otsu*/
int cout_otsu(Matrix Hist, int i) {
  int ** hist = MatGetInt(Hist);
  int s = 0;
  int mu1 = 0;
  int mu2 = 0;
  int q1 = 0;
  int q2 = 0;
  for(int k = 0; k <= 255; k++) {
    if(k < i) {
      q1 += hist[0][k];
      mu1 += k * hist[0][k];
    } else {
      q2 += hist[0][k];
      mu2 += k * hist[0][k];
    }
  }
  
  mu1 = q1 ? mu1/q1 : 0;
  mu2 = q2 ? mu2/q2 : 0;

  for(int k = 0; k <= 255; k++) {
    if (k <= i-1) {
      s += hist[0][k]*(k-mu1)*(k-mu1);
    } else {
      s += hist[0][k]*(k-mu2)*(k-mu2);
    }
  }
  return s;
}
/*---------------------------------*/

/*
 * Calcul du seuil d'Otsu
 * Entrée : histogramme (matrice de int 1 x 256)
 * Sortie : seuil (niveau de gris) obtenu par la méthode d'Otsu
 */
unsigned char Otsu(Matrix Hist)
{
  int S = 1;
  int min = cout_otsu(Hist, 1);
  for(int i = 2; i < 255; i++) {
    int cout_i = cout_otsu(Hist,i);
    if (cout_i < min) {
      min = cout_i;
      S = i;
    }
    //printf("%f\n",min);
  }
  return S;
}




/*
 * Calcul de l'histogramme cumulé à partir de l'histogramme
 * Entrée : histogramme (matrice de int 1 x 256)
 * Sortie : histogramme cumulé (matrice de int 1 x 256)
 */
Matrix Hist2CumHist(Matrix Hist)
{
  AFAIRE(NULL);
}


/*
 * Application d'une transformation ponctuelle à une image de niveaux de gris
 * Entreés : image initiale (niveaux de gris) et
 * transformation ponctuelle (matrice de int 1 x 256)
 * Sortie : image de niveaux de gris transformée
 */
Image AppLUT(Image Im, Matrix LUT)
{
  AFAIRE(NULL);
}


/*
 * Spécification d'histogramme
 * Entrées : histogramme cumulé de l'image et histogramme cumulé desiré
 * (on suppose que le dernier élément des deux histogrammes cumulés sont
 * les mêmes, c'est-à-dire qu'ils décrivent des images contenant le même nombre
 * de pixels)
 * Sortie : transformation ponctuelle (matrice 1 x 256)
 */
Matrix HistSpecif(Matrix CumHist, Matrix DesCumHist)
{
  AFAIRE(NULL);
}


/*
 * Vérification de la validité d'une matrice représentant un élément
 * structurant binaire (pour l'érosion, la dilatation, etc.)
 * Entrée : matrice représentant un élément structurant
 * Sortie : 0 si la matrice est valide,
            SE_NOT_ODD si son nombre de lignes ou de colonnes n'est pas impair
            SE_NOT_INT si elle ne contient pas que des entiers
            SE_NOT_BIN si elle ne contient pas que des 0 et des 1
*/
int NotValidBinSE(Matrix StructuringElement)
{
  int **ES,NbLig,NbCol,i,j;

  if (MatType(StructuringElement)!=Int)
    return SE_NOT_INT;
  NbLig=MatNbRow(StructuringElement);
	if ((NbLig%2)!=1)
	  return SE_NOT_ODD;
  NbCol=MatNbCol(StructuringElement);
	if ((NbCol%2)!=1)
	  return SE_NOT_ODD;
  ES=MatGetInt(StructuringElement);
  for (i=0;i<NbLig;i++)
    for (j=0;j<NbCol;j++)
      if (ES[i][j]!=0 && ES[i][j]!=1)
        return SE_NOT_BIN;
  return 0;
}


/*
 * Vérification de la validité d'une matrice représentant un élément
 * structurant ternaire (pour la transformation "tout ou rien")
 * Entrée : matrice représentant un élément structurant
 * Sortie : 0 si la matrice est valide,
            SE_NOT_ODD si son nombre de lignes ou de colonnes n'est pas impair
            SE_NOT_INT si elle ne contient pas que des entiers
            SE_NOT_TERN si elle ne contient pas que des 0, des 1 et des 2
*/
int NotValidTernSE(Matrix StructuringElement)
{
  int **ES,NbLig,NbCol,i,j;

  if (MatType(StructuringElement)!=Int)
    return SE_NOT_INT;
  NbLig=MatNbRow(StructuringElement);
	if ((NbLig%2)!=1)
	  return SE_NOT_ODD;
  NbCol=MatNbCol(StructuringElement);
	if ((NbCol%2)!=1)
	  return SE_NOT_ODD;
  ES=MatGetInt(StructuringElement);
  for (i=0;i<NbLig;i++)
    for (j=0;j<NbCol;j++)
      if (ES[i][j]!=0 && ES[i][j]!=1 && ES[i][j]!=2)
        return SE_NOT_TERN;
  return 0;
}


/*
 * Amincissement d'une image binaire
 * Entreés : image binaire initiale et élément structurant (matrice de int
 * contenant uniquement des 0, des 1 et des 2 signifiant "peu importe")
 * Sortie : image binaire transformée
 */
Image Thinning(Image Im, Matrix StructuringElement)
{
  AFAIRE(NULL);
}

