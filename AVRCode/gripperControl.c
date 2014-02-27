/*
 * GccApplication1.c
 *
 * Created: 14.02.2014 15:28:02
 *  Author: vegardvo
 */ 


#include <avr/io.h>

#define PWMPIN PB0 // PWM out on pin 0	
#define SWT PB1	   // Button input on pin 1

int main(void)
{
	// Set pin directions
	DDRB &= ~(1 << SWT);
	DDRB |= (1 << PWMPIN);
	
	// Enable PWM
	TCCR0A |= (1 << COM0A1)  | (1 << WGM01) | (1 << WGM00);  // Fast PWM
	TCCR0B |= (0 << WGM02);									 // Fast PWM
						
	
	// Set prescaler
	TCCR0B |= (1 << CS01) | (1 << CS00); // Clock select 64 prescaler
	
	// Init PWM signal
	OCR0A = 0;
    while(1)
    {

//		OCR0A = 16;
		OCR0A = 31;
		//OCR0A = 0;
}