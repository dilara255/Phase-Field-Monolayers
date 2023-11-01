#pragma once

#include "miscStdHeaders.h"
#include "miscDefines.hpp"
#include "core.hpp"

//TODO: template
//TODO: why not use std::bitset instead?
//TODO: if not, rethink counting method: 
//-- https://en.cppreference.com/w/cpp/numeric/popcount
//--- But Byte-At-A-Time Lookup Table seems to be within ~17% of native popcount
//---- ( https://github.com/intvnut/bit_count )

namespace AZ {

	class FlagField32 {
	
	public:
		FlagField32() { }

		void loadField(uint32_t flagsToLoad) {
			flags = flagsToLoad;
			updateHowManyAreOn();
		}

		uint32_t getField() const {
			return flags;
		}

		void updateHowManyAreOn() {
			areOn = 0;

			for (int i = 0; i < (int)(sizeof(flags) * BITS_IN_A_BYTE); i++) {
					areOn += ((flags & BIT(i)) != 0);
			}
		}
		
		int howManyAreOn() const {
			return areOn;
		}

		bool isBitOn(uint32_t bit) {
			if (bit > 31) return false; //doesn't exist: not active

			return (flags & BIT(bit));
		}

		void setBitOn(int bit) {
			if (bit > 31) return; //doesn't exist

			areOn += ((flags & BIT(bit)) == 0);

			flags = flags | (BIT(bit));
		}

		void setBitOff(int bit) {
			if (bit > 31) return; //doesn't exist

			areOn -= ((flags & BIT(bit)) > 0);

			flags = flags & (~(BIT(bit)));
		}

		std::string getFlagsAsString() {
			std::bitset bits = std::bitset<4*BITS_IN_A_BYTE>(flags);
			return bits.to_string();
		}

		/*TODO: implement toggle
		bool toggleBit(int bit) {

		}
		*/

	private:
		uint32_t flags = 0;
		int areOn = 0;
	};

	class FlagField128 {

	public:
		FlagField128() { flags[0] = 0; flags[1] = 0; flags[2] = 0; flags[3] = 0; }

		bool loadField(uint32_t flagsToLoad, uint32_t toWhichField) {
			if (toWhichField > 3) return false; //doesn't exist: not active

			flags[toWhichField] = flagsToLoad;

			updateHowManyAreOn();

			return true;
		}

		uint32_t getField(uint32_t whichField) const {
			if (whichField > 3) return 0; //doesn't exist: not active

			return flags[whichField];
		}

		void updateHowManyAreOn() {
			areOn = 0;
			for (int i = 0; i < blocks; i++) {
				for (int j = 0; j < (int)(sizeof(flags[0]) * BITS_IN_A_BYTE); j++) {
					areOn += ( (flags[i] & BIT(j)) != 0);
				}
			}
		}

		int howManyAreOn() const {
			return areOn;
		}

		bool isBitOn(int bit) {
			if (bit > 127) return false; //doesn't exist: not active

			int block = bit / (sizeof(flags[0]) * BITS_IN_A_BYTE);
			int bitOnTheBlock = bit - (block * (sizeof(flags[0]) * BITS_IN_A_BYTE));

			return (flags[block] & BIT(bitOnTheBlock));
		}

		//TODO: maybe extract calculation og block and bitOnTheBlock?
		void setBitOn(int bit) {
			if (bit > 127) return; //doesn't exist

			int block = bit / (sizeof(flags[0]) * BITS_IN_A_BYTE);
			int bitOnTheBlock = bit - (block * (sizeof(flags[0]) * BITS_IN_A_BYTE));

			//if bit was off:
			areOn += ((flags[block] & BIT(bitOnTheBlock)) == 0);

			flags[block] = flags[block] | (BIT(bitOnTheBlock));
		}

		void setBitOff(int bit) {
			if (bit > 127) return; //doesn't exist
			
			int block = bit / (sizeof(flags[0]) * BITS_IN_A_BYTE);
			int bitOnTheBlock = bit - (block * (sizeof(flags[0]) * BITS_IN_A_BYTE));

			//if bit was on:
			areOn -= ((flags[block] & BIT(bitOnTheBlock)) > 0);

			flags[block] = flags[block] & (~(BIT(bitOnTheBlock)));
		}

		int getNumberOfBlocks() {
			return blocks;
		}

		size_t sizeOfBlockInBytes() {
			return sizeof(flags[0]);
		}

		std::string getFlagsAsString(int whichField) {
			std::bitset bits = std::bitset<4*BITS_IN_A_BYTE>(flags[whichField]);
			return bits.to_string();
		}

		/*TODO: implement toggle
		void toggleBit(int bit) {

		}
		*/

	private:
		int blocks = 4;
		uint32_t flags[4];
		int areOn = 0;
	};
}
