public class UTMConverter {

	final static int a = 6378137;
	final static double f = 1.0 / 298.257223563;
	final static double n = f / (2 - f);
	final static double k0 = 0.9996;
	final static double E0 = 500e3;
	final static double N0 = 10000e3;
	final static double e = Math.sqrt(f * (2 - f));
	final static double n2 = n * n;
	final static double n3 = n * n2;
	final static double n4 = n * n3;
	final static double n5 = n * n4;
	final static double n6 = n * n5;

	public static double atanh(double x) {
		return 0.5 * (Math.log(1 + x) - Math.log(1 - x));
	}

	public static double asinh(double x) {
		return Math.log(x + Math.sqrt(1 + x * x));
	}

	public static double[] converter(double lat, double lon) {

		int zone = (int) Math.floor((lon + 180) / 6) + 1;
		double λ0 = Math.toRadians((zone - 1) * 6 - 180 + 3);

		double φ = Math.toRadians(lat);
		double λ = Math.toRadians(lon) - λ0;

		double A = a / (1 + n) * (1.0 + n2 / 4 + n4 / 64 + n6 / 256);

		double[] α = new double[] { 0, 1.0 / 2 * n - 2.0 / 3 * n2 + 5.0 / 16 * n3,
				13.0 / 48 * n2 - 3.0 / 5 * n3, 61.0 / 240 * n3 };

		double sinφ = Math.sin(φ);
		double sinλ = Math.sin(λ);
		double cosλ = Math.cos(λ);

		double nt = 2.0 * Math.sqrt(n) / (1 + n);
		double t = Math.sinh(atanh(sinφ) - nt * atanh(nt * sinφ));
		double ξʹ = Math.atan(t / cosλ);
		double ηʹ = atanh(sinλ / Math.sqrt(1 + t * t));
		//double τ = Math.tan(φ);
		//double σ = Math.sinh(e * atanh(e * τ / Math.sqrt(1 + τ * τ)));
		//double τʹ = τ * Math.sqrt(1 + σ * σ) - σ * Math.sqrt(1 + τ * τ);
		//double ξʹ = Math.atan2(τʹ, cosλ);
		//double ηʹ = asinh(sinλ / Math.sqrt(τʹ * τʹ + cosλ * cosλ));

		double η = ηʹ;
		for (int j = 1; j <= 3; j++)
			η += α[j] * Math.cos(2.0 * j * ξʹ) * Math.sinh(2.0 * j * ηʹ);

		double ξ = ξʹ;
		for (int j = 1; j <= 3; j++)
			ξ += α[j] * Math.sin(2.0 * j * ξʹ) * Math.cosh(2.0 * j * ηʹ);

		double x = k0 * A * η;
		double y = k0 * A * ξ;

		x = E0 + x;
		if (y < 0)
			y = N0 + y;

		System.out.println("x=" + x);
		System.out.println("y=" + y);
		return (new double[] { x, y });
	}

	public static void main(String[] args) {

		converter(48.8582, 2.2945);
		converter(32.744353, -117.249111);
		System.out.println("aaa");

	}

}
