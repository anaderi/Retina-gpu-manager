#include "Tools.h"

int* h_no_sensors;
int* h_no_hits;
int* h_sensor_Zs;
int* h_sensor_hitStarts;
int* h_sensor_hitNums;
int* h_hit_IDs;
float* h_hit_Xs;
float* h_hit_Ys;
int* h_hit_Zs;


// struct PixelEvent {
//   int noSensors;
//   int noHits;
//   std::vector<int> sensorZs;
//   std::vector<int> sensorHitStarts;
//   std::vector<int> sensorHitsNums;
//   std::vector<int> hitIDs;
//   std::vector<float> hitXs;
//   std::vector<float> hitYs;
//   std::vector<int> hitZs;
// };

void setHPointersFromPixelEvent(const PixelEvent& event){
	h_no_sensors = (int*) &event.noSensors;
	h_no_hits = (int*) &event.noHits;
	h_sensor_Zs = (int*) &event.sensorZs[0];
	h_sensor_hitStarts = (int*) &event.sensorHitStarts[0];
	h_sensor_hitNums = (int*) &event.sensorHitsNums[0];
	h_hit_IDs = (int*) &event.hitIDs[0];
	h_hit_Xs = (float*) &event.hitXs[0];
	h_hit_Ys = (float*) &event.hitYs[0];
	h_hit_Zs = (int*) &event.hitZs[0];
}

void setHPointersFromInput(char*& input){
	h_no_sensors = (int*) &input[0];
	h_no_hits = (int*) (h_no_sensors + 1);
	h_sensor_Zs = (int*) (h_no_hits + 1);
	h_sensor_hitStarts = (int*) (h_sensor_Zs + h_no_sensors[0]);
	h_sensor_hitNums = (int*) (h_sensor_hitStarts + h_no_sensors[0]);
	h_hit_IDs = (int*) (h_sensor_hitNums + h_no_sensors[0]);
	h_hit_Xs = (float*) (h_hit_IDs + h_no_hits[0]);
	h_hit_Ys = (float*) (h_hit_Xs + h_no_hits[0]);
	h_hit_Zs = (int*) (h_hit_Ys + h_no_hits[0]);
}

void printInfo(){
	std::cout << "Read info:" << std::endl
		<< " no sensors: " << h_no_sensors[0] << std::endl
		<< " no hits: " << h_no_hits[0] << std::endl
		<< "First 5 sensors: " << std::endl;

	for (int i=0; i<5; ++i){
		std::cout << " Zs: " << h_sensor_Zs[i] << std::endl
			<< " hitStarts: " << h_sensor_hitStarts[i] << std::endl
			<< " hitNums: " << h_sensor_hitNums[i] << std::endl << std::endl;
	}

	std::cout << "First 5 hits: " << std::endl;

	for (int i=0; i<5; ++i){
		std::cout << " hit_id: " << h_hit_IDs[i] << std::endl
			<< " hit_X: " << h_hit_Xs[i] << std::endl
			<< " hit_Y: " << h_hit_Ys[i] << std::endl
			<< " hit_Z: " << h_hit_Zs[i] << std::endl << std::endl;
	}
}

void readFile(std::string filename, char*& input, int& size){
	// Give me them datas!!11!
	std::ifstream infile (filename.c_str(), std::ifstream::binary);

	// get size of file
	infile.seekg(0, std::ifstream::end);
	size = infile.tellg();
	infile.seekg(0);

	// read content of infile with pointers
	input = (char*) malloc(size);
	infile.read (input, size);
	infile.close();

	setHPointersFromInput(input);
}

void quickSortInput(char*& input){
	for(int i=0; i<h_no_sensors[0]; i++)
        quickSort(h_hit_Xs, h_hit_Ys, h_hit_IDs, h_hit_Zs,
		    h_sensor_hitStarts[i], h_sensor_hitStarts[i] + h_sensor_hitNums[i]);
}

void quickSort(float*& hit_Xs, float*& hit_Ys, int*& hit_IDs, int*& hit_Zs, int _beginning, int _end)
{
	const int max_levels = 300;
	int beg[max_levels], end[max_levels], i=0, L, R, swap;

	double piv, d1;
	int i1, i2;

	beg[0]=_beginning; end[0]=_end;
	while (i>=0) {
		L=beg[i]; R=end[i]-1;
		if (L<R) {

			piv = hit_Xs[L];
			d1  = hit_Ys[L];
			i1  = hit_IDs[L];
			i2  = hit_Zs[L];

			while (L<R) {
				while (hit_Xs[R] >= piv && L < R) R--;
				if (L<R){
					hit_Xs[L] = hit_Xs[R];
					hit_Ys[L] = hit_Ys[R];
					hit_Zs[L] = hit_Zs[R];
					hit_IDs[L] = hit_IDs[R];
					L++;
				}

				while (hit_Xs[L] <= piv && L < R) L++;
				if (L<R){
					hit_Xs[R] = hit_Xs[L];
					hit_Ys[R] = hit_Ys[L];
					hit_Zs[R] = hit_Zs[L];
					hit_IDs[R] = hit_IDs[L];
					R--;
				}
			}
			hit_Xs[L] = piv;
			hit_Ys[L] = d1;
			hit_IDs[L] = i1;
			hit_Zs[L] = i2;

			beg[i+1]=L+1; end[i+1]=end[i]; end[i++]=L;
			if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
				swap=beg[i]; beg[i]=beg[i-1]; beg[i-1]=swap;
				swap=end[i]; end[i]=end[i-1]; end[i-1]=swap;
			}
		}
		else {
			i--;
		}
	}
}
