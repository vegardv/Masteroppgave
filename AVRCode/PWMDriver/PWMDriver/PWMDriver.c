/*
 * PWMDriver.c
 *
 * Created: 07.04.2014 13:08:07
 *  Author: vegardvo
 */ 


#include <avr/io.h>

#define PWMPIN PB0 // PWM out on pin 0
#define EN PB1 // Enable input on pin 1
#define POS PB2	   // Position input on pin 2

int main(void)
{
		// Set pin directions
		DDRB &= ~(1 << EN);
		DDRB &= ~(1 << POS);
		//DDRB |= (1 << PWMPIN);
		
		// Enable pull-up
		//PORTB |= (1 << EN);
		PORTB |= (1 << POS);
		//MCUCR |= (1 << PUD);
		
		// Enable PWM
		TCCR0A |= (1 << COM0A1)  | (1 << WGM01) | (1 << WGM00);  // Fast PWM
		TCCR0B |= (0 << WGM02);									 // Fast PWM
		
		
		// Set prescaler
		TCCR0B |= (1 << CS01) | (1 << CS00); // Clock select 64 prescaler
		
		// Init PWM signal
		OCR0A = 0;
    while(1)
    {
		// Read input
		// Enabled high
		// Pos high -> close
	    if(PINB & 0b00000010)
		{
			DDRB |= (1 << PWMPIN);
			if(!(PINB & 0b00000100))
			{
				OCR0A = 16;
			}
			else
			{
				OCR0A = 31;
			}
						
		}
		else
		{
			//OCR0A = 0;
			DDRB &= ~(1 << PWMPIN);
	    }
	}
		
}