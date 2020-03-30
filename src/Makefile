CFLAGS = -g
CFLAGS = -O3 -march=native
LDLIBS += -lm -lhdf5_serial

OBJ = iio.o
BIN = splat h5ls h5get localize matchl12 mscalar h5v h5getinfo

all : $(BIN)

$(BIN) : $(OBJ)

clean: ; $(RM) -f $(BIN) $(OBJ)
