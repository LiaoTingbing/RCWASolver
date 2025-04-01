#include "common.h"


class Device
{
public:
	Device();
	~Device();

private:
	cx_mat* Index;
	cx_mat ER;
	colvec Lnode;


};

Device::Device()
{
}

Device::~Device()
{
}