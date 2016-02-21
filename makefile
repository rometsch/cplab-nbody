################################################################################
# Automatically-generated file. Do not edit!
################################################################################

RM := rm -rf

# All of the sources participating in the build are defined here
CPP_SRCS := \
./gnuplot_i.hpp \
./Integrators.cpp \
./Logger.cpp \
./main.cpp \
./nbsys.cpp \

OBJS += \
./Integrators.o \
./Logger.o \
./main.o \
./nbsys.o \

CPP_DEPS += \
./Integrators.d \
./Logger.d \
./main.d \
./nbsys.d \

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: n-body

# Tool invocations
n-body: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: MacOS X C++ Linker'
	g++ -O3 -o "n-body" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '
	mkdir aufgabe2a
	mkdir aufgabe2b
	mkdir aufgabe3

# Other Targets
clean:
	-$(RM) $(C++_DEPS)$(OBJS)$(C_DEPS)$(CC_DEPS)$(CPP_DEPS)$(EXECUTABLES)$(CXX_DEPS)$(C_UPPER_DEPS) n-body
	-@echo ' '
