using System;
using System.IO;
using System.Diagnostics;

namespace homework
{

	public class node<T>
	{
		public T item;
		public node<T> next;
		public node(T item) { this.item = item; }

	}

	public class list<T>
	{
		public node<T> first = null, current = null;
		public list() {}
		public void push(T item)
		{
			if (first == null)
			{
				first = new node<T>(item);
				start();
			}
			else
			{
				while (current.next != null) next();
				current.next = new node<T>(item);
				next();
			}
		}
		public void start() { current = first; }
		public void next() { current = current.next; }

		public node<T> get_node(int i)
		{
			start();
			for (int j = 0; j < i && current != null; j++) next();
			return current;
		}

		public T get(int i)
		{
			node<T> res = get_node(i);
			Debug.Assert(res != null, "Index out of bounds");
			return res.item;

		}

		public void remove(int i)
		{
			if (i == 0) first = first.next; 
			else
            {
				node<T> itnode = get_node(i - 1);
				itnode.next = itnode.next.next;
			}
		}

		public int Length
        {
			get
			{
				start();
				int res = 0;
				while (current.next != null) { next(); res++; }
				return res;
			}
        }
	}

	public class list_methods
	{

		public static list<vector> file_to_list(string arg)
		{
			string workingDirectory = Environment.CurrentDirectory;
			string location = Directory.GetParent(workingDirectory).Parent.Parent.FullName + "\\" + arg;
			var list = new list<genlist<double>>();
			var list_res = new list<vector>();
			int listlen = 0;
			int lineno = 0;
			using (StreamReader w = new StreamReader(location))
			{
				char[] delimiters = { ' ', '\t', ','};
				var options = StringSplitOptions.RemoveEmptyEntries;
				
				for (string line = w.ReadLine(); line != null; line = w.ReadLine())
				{
					string[] words = line.Split(delimiters, options);
					int n = words.Length;
					double[] nums = new double[n];
					bool success = true;
					for (int i = 0; i < n; i++)
                    {
						success = double.TryParse(words[i], out nums[i]);
						if (!success) break;
					}
					if (success)
                    {
						while (listlen < n)
						{
							genlist<double> lst = new genlist<double>();
							for (int i = 0; i < lineno; i++) lst.push(0);
							list.push(lst);
							listlen++;
						}

						for (int i = 0; i < n; i++) list.get(i).push(nums[i]);
						for (int i = n; i < listlen; i++) list.get(i).push(0);
						lineno++;
					}
				}
			}
			for (int i = 0; i < listlen; i++)
			{
				vector v = new vector(lineno);
				v.setData(list.get(i).ToArray());
				list_res.push(v);
			}

			return list_res;
		}
		 
		public static void test()
		{
			list<int> a = new list<int>();
			a.push(1);
			a.push(2);
			a.push(3);
			a.remove(1);
			for (a.start(); a.current != null; a.next())
			{
				Debug.WriteLine(a.current.item);
			}
		}
	}

	public class genlist<T>
	{
		public T[] data;
		public int capacity => data.Length; // property
		public int size;
		public genlist(int capacity = 0) {
			data = new T[capacity];
			this.size = 0;
		}

		public T this[int i]
		{ 
			get // indexer
			{
				Debug.Assert(i < size, "Index out of bounds");
				return data[i];
			}
			set // setter: for example, A[i]=5;
			{
				Debug.Assert(i < size, "Index out of bounds");
				data[i] = value;
			}
		}

		public void push(T item)
		{ /* add item to list */
			if (capacity == 0) data = new T[1];
			else if (size == capacity)
			{
				T[] newdata = new T[2 * capacity];
				for (int i = 0; i < size; i++) newdata[i] = data[i];
				data = newdata;
			}
			data[size] = item;
			size++;
			
			
		}

		public void push_series(T[] items)
		{
			foreach (T item in items) push(item);
		}

		public void remove(int i)
        {
			for (int j = i; j < size - 1; j++) data[j] = data[j + 1];
			size--;
		}

		public T[] ToArray()
        {
			T[] arr = new T[size];
			for (int i = 0; i < size; i++) arr[i] = data[i];
			return arr;
		}
	}
}
